#!/usr/bin/env Rscript

# Arguments from bash
print(commandArgs(trailingOnly = TRUE))

args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) == 2)){
    stop("usage: refine_passes.R config.yaml scan2_object_update.rda")
}



## Remove before calling from bash ##
#yaml_fn = '/home/ch253794/workdir/PROJECTS/SCAN2_FAST/REFINE/scan.yaml'
#scan2_fn = '/lab-share/Gene-Walsh-e2/Public/indData/emre/PROJECTS/NHP_MUTATION/SCAN2/CHIMP/PTA_6237/call_mutations/6237-C8-Exc/scan2_object.rda'
#samp_bulk_bam = '/home/ch253794/nr_workdir/PROJECTS/NHP_MUTATION/BULK/CHIMP_6237/VARIANTS/6237_base_recalib.bam'
#metafile = '/home/ch253794/workdir/PROJECTS/NHP_MUTATION/SCAN2/METADATA/panel_chimp.csv'
#id_frac = 0.5
#gb = 'panTro6'
#bsgenome = 'BSgenome.Ptroglodytes.UCSC.panTro6'
#fastafile = '/home/ch253794/workdir/REFERENCES/BWA_REFERENCES/panTro6/panTro6.fa'

# Get from bash
yaml_fn = args[1]
scan2_fn = args[2]

# Extract from config file
config = yaml::read_yaml(yaml_fn)
samp_bulk_bam = config$bulk_bam
gb = config$genome
bsgenome = config$bsgenome
fastafile = config$ref
id_frac = config$indel_fraction_delta
metafile = config$makepanel_metadata
allpasses_fn = config$allpasses
allartifs_fn = config$allartifs


if(gb == 'hs37d5'){
	gb = 'hg19'
	print('Changing hs37d5 to hg19; it does not change the behavior for this function')
}

# Load required libraries
suppressPackageStartupMessages(suppressWarnings(require(scanZooTrial1)))
suppressPackageStartupMessages(suppressWarnings(require(dplyr)))
suppressPackageStartupMessages(suppressWarnings(require(ggpubr)))
suppressPackageStartupMessages(suppressWarnings(require(parallel)))
suppressPackageStartupMessages(suppressWarnings(require(lme4)))
suppressPackageStartupMessages(suppressWarnings(require(cowplot)))
suppressPackageStartupMessages(suppressWarnings(require(MutationalPatterns)))
suppressPackageStartupMessages(suppressWarnings(require(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(require(Seurat)))
suppressPackageStartupMessages(suppressWarnings(require(vcfR)))
suppressPackageStartupMessages(suppressWarnings(require(data.table)))

# For now, source the scripts that extract info from scan2 objects
scripts = paste0(config$scripts, '/UTILITY_FUNCTIONS.R')
source(scripts)

####
## PREPARE FOR RUN
####

# Derived variables
allpasses = readRDS(allpasses_fn)
pass_artif = readRDS(allartifs_fn)
meta = read.csv(metafile)

outdir = 'objects'
dir.create(outdir, showWarnings=F)

tmpdir = 'tmpdir'
dir.create(tmpdir, showWarnings=F)

# Adjustments
allpasses$id = paste0(allpasses$chr, '_', allpasses$pos, '_', allpasses$refnt, '_', allpasses$altnt)

# Load gatk
res = get(load(scan2_fn))
gatk = res@gatk
gatk$id = paste(gatk$chr, gatk$pos, gatk$refnt, gatk$altnt, sep = '_')
gatk_indel = gatk[gatk$muttype == 'indel',]
gatk_snv = gatk[gatk$muttype == 'snv',]

rm(gatk)
gc()

# Load all passes and artifacts
pass_artif_snv = pass_artif[pass_artif$muttype == 'snv']
pass_artif_indel = pass_artif[pass_artif$muttype == 'indel']
allpasses_snv = allpasses[allpasses$muttype == 'snv']
allpasses_indel = allpasses[allpasses$muttype == 'indel']

rm(pass_artif, allpasses)
gc()

####
## INDELS
####

## Remove if in bulk ##
samp = basename(dirname(scan2_fn))
print(paste0('Processing sample for indel ', samp))

## Add the ones that are in the bulk up to 3 reads ##
samp_passes_indel1 = gatk_indel[cigar.id.test & cigar.hs.test & dp.test &
	abc.test & min.sc.alt.test & csf.test & bulk.gt == '0/0' &
	mda.fdr < 0.01 & lysis.fdr < 0.01 & csf.test & balt.lowmq <= 2,]
print(paste0('Total indel passes ', sum(gatk_indel$pass)))
print(paste0('Total indel passes with relaxed balt.lowmq ', nrow(samp_passes_indel1)))

## GATK - Regenotype new passes, remove if present in bulk ##
sites = paste0(samp_passes_indel1$chr, ':', samp_passes_indel1$pos, '-', samp_passes_indel1$pos)
outvcf = paste0(outdir, '/Indel_Regenotype_', samp, '.vcf')

gatk_args <- c(
  "HaplotypeCaller",
  "-R", fastafile,
  "-I", samp_bulk_bam,
  "-O", outvcf,
  "--minimum-mapping-quality", "0",
  "--max-reads-per-alignment-start", "0",
  "--dont-use-soft-clipped-bases", "true",
  "--min-base-quality-score", "0",
  "--emit-ref-confidence", "GVCF",
  "--standard-min-confidence-threshold-for-calling", "0",
  "--pcr-indel-model", "NONE",
  "--min-pruning", "0",
  "--interval-padding", "100",
  "--min-dangling-branch-length", "0",
  "--allow-non-unique-kmers-in-ref", "true",
  unlist(lapply(sites, function(x) c("-L", x)))
)

# Run the command
system2("gatk", args = gatk_args)
indel_regen = read.vcfR(outvcf)
indel_regen_fix = indel_regen@fix

# Check if alt allele is present
alts = strsplit(indel_regen_fix[, 'ALT'], ',')

regen_ids = lapply(1:length(alts), function(k){paste(indel_regen_fix[k, 'CHROM'], indel_regen_fix[k, 'POS'], indel_regen_fix[k, 'REF'], alts[[k]], sep = '_')}) %>% unlist

filt_bulkgatk_indel = samp_passes_indel1[samp_passes_indel1$id %in% regen_ids,]

## Remove calls in sites with many indels and clippings ##
gatk_indel$sc_id_frac = (gatk_indel$ID.cigars/gatk_indel$dp.cigars) - (gatk_indel$ID.cigars.bulk/gatk_indel$dp.cigars.bulk)
gatk_indel$sc_id_frac_afNorm = gatk_indel$sc_id_frac / gatk_indel$af
indel_cigarpass = gatk_indel[sc_id_frac_afNorm > id_frac, id]

filt_morebulk_indel = samp_passes_indel1[!(samp_passes_indel1$id %in% indel_cigarpass),]


## Remove if multiple indels are nearby ##
filt_cls_indel = find_clustered_snvs(obj = res, inpass = samp_passes_indel1)


## Remove if recurrent, even if filtered as a pre-amplification artifact ##
# Check if any two cells share unusually high number of mutations
subj = meta[meta$sample == samp, 'donor']
freq = allpasses_indel[allpasses_indel$donor == subj, 'id'] %>% table
freq_two = names(freq[freq == 2])
freq_two_cnts = allpasses_indel[allpasses_indel$id %in% freq_two, 'Sample'] %>%
			table %>% sort(., decreasing = T)
shared_ratio = (min(freq_two_cnts[1:2]) - freq_two_cnts[3]) / sum(na.omit(gatk_indel$pass))
shared_ratio = ifelse(is.na(shared_ratio), 0, shared_ratio)

# If very few sSNVs are recurring, filter any recurrent sSNVs
# If many sSNVs are recurring between two samples, only filter other recurrent sSNVs
if(shared_ratio < 0.2){
	pass_artif_other = pass_artif_indel[pass_artif_indel$Sample != samp & pass_artif_indel$pass == F, ]
	pass_artif_other$chr = gsub('chr', '', pass_artif_other$chr)
	recurrent_artif = names(which(table(pass_artif_other$id) >= 1))
} else{
	other_sample = setdiff(names(freq_two_cnts[1:2]), samp)
	pass_artif_other = pass_artif_indel[!(pass_artif_indel$Sample %in% c(samp, other_sample))
					& pass_artif_indel$pass == F, ]
	pass_artif_other$chr = gsub('chr', '', pass_artif_other$chr)
	recurrent_artif1 = names(which(table(pass_artif_other$id) >= 1))

	pass_artif_other = pass_artif_indel[pass_artif_indel$Sample != samp &
				pass_artif_indel$pass == F, ]
	pass_artif_other$chr = gsub('chr', '', pass_artif_other$chr)
	recurrent_artif2 = names(which(table(pass_artif_other$id) >= 2))

	recurrent_artif = union(recurrent_artif1, recurrent_artif2)
}

filt_recur_indel = samp_passes_indel1[samp_passes_indel1$id %in% recurrent_artif,]

# Assign passes
filters_indel = do.call(rbind, list(filt_bulkgatk_indel, filt_morebulk_indel, filt_cls_indel, filt_recur_indel))
samp_passes_indel2 = samp_passes_indel1[!(samp_passes_indel1$id %in% filters_indel$id),]
gatk_indel$pass = F
gatk_indel[id %in% samp_passes_indel2$id, pass := T]

####
## SNVs
####

print(paste0('Processing sample for sSNV ', samp))

## Add the ones that are in the bulk up to 2 reads ##
samp_passes_snv1 = gatk_snv[cigar.id.test & cigar.hs.test & dp.test &
	abc.test & min.sc.alt.test & csf.test & bulk.gt == '0/0' &
	mda.fdr < 0.01 & lysis.fdr < 0.01 & csf.test & balt.lowmq <= 2,]

samp_passes_snv1 = add_mutsigs(mutdf = samp_passes_snv1,
			bsgenome=bsgenome,
			gb=gb)

if(gb == 'hg19'){
	samp_passes_snv1$chr = gsub('chr', '', samp_passes_snv1$chr)
}

print(paste0('Total snv passes ', sum(gatk_snv$pass)))
print(paste0('Total snv passes with relaxed balt.lowmq ', nrow(samp_passes_snv1)))

## GATK - Regenotype new passes, remove if present in bulk ##
sites = paste0(samp_passes_snv1$chr, ':', samp_passes_snv1$pos, '-', samp_passes_snv1$pos)
outvcf = paste0(outdir, '/Snv_Regenotype_', samp, '.vcf')

gatk_args <- c(
  "HaplotypeCaller",
  "-R", fastafile,
  "-I", samp_bulk_bam,
  "-O", outvcf,
  "--minimum-mapping-quality", "0",
  "--max-reads-per-alignment-start", "0",
  "--dont-use-soft-clipped-bases", "true",
  "--min-base-quality-score", "0",
  "--emit-ref-confidence", "GVCF",
  "--standard-min-confidence-threshold-for-calling", "0",
  "--min-pruning", "0",
  "--interval-padding", "10",
  "--min-dangling-branch-length", "0",
  "--allow-non-unique-kmers-in-ref", "true",
  unlist(lapply(sites, function(x) c("-L", x)))
)

# Run the command
system2("gatk", args = gatk_args)
snv_regen = read.vcfR(outvcf)
snv_regen_fix = snv_regen@fix

# Check if alt allele is present
alts = strsplit(snv_regen_fix[, 'ALT'], ',')
regen_ids = lapply(1:length(alts), function(k){paste(snv_regen_fix[k, 'CHROM'], snv_regen_fix[k, 'POS'], snv_regen_fix[k, 'REF'], alts[[k]], sep = '_')}) %>% unlist
filt_bulkgatk_snv = samp_passes_snv1[samp_passes_snv1$id %in% regen_ids,]


## MPILEUP - Regenotype new passes, remove if present in bulk ##
samp_passes_snv1 = add_mutsigs(mutdf = samp_passes_snv1,
			bsgenome=bsgenome, gb=gb)

if(gb == 'hg19'){
	samp_passes_snv1$chr = gsub('chr', '', samp_passes_snv1$chr)
}

snv_bulk_passes = remove_bulk_mpileup(passes_file = samp_passes_snv1,
			bulk_bam_file = samp_bulk_bam,
			outfile = tempfile(tmpdir = tmpdir, fileext = ".txt"))
filt_bulkmp_snv = samp_passes_snv1[!(samp_passes_snv1$id %in% snv_bulk_passes$id),]

## Remove if multiple sSNVs are nearby ##
filt_cls_snv = find_clustered_snvs(obj = res, inpass = samp_passes_snv1)

## Remove if recurrent, even if filtered as a pre-amplification artifact ##
# Check if any two cells share unusually high number of mutations
subj = meta[meta$sample == samp, 'donor']
freq = allpasses_snv[allpasses_snv$donor == subj, 'id'] %>% table
freq_two = names(freq[freq == 2])
freq_two_cnts = allpasses_snv[allpasses_snv$id %in% freq_two, 'Sample'] %>%
			table %>% sort(., decreasing = T)
shared_ratio = (min(freq_two_cnts[1:2]) - freq_two_cnts[3]) / sum(na.omit(gatk_snv$pass))
shared_ratio = ifelse(is.na(shared_ratio), 0, shared_ratio)

# If very few sSNVs are recurring, filter any recurrent sSNVs
# If many sSNVs are recurring between two samples, only filter other recurrent sSNVs
if(shared_ratio < 0.2){
	pass_artif_other = pass_artif_snv[pass_artif_snv$Sample != samp & pass_artif_snv$pass == F, ]
	pass_artif_other$chr = gsub('chr', '', pass_artif_other$chr)
	recurrent_artif = names(which(table(pass_artif_other$id) >= 1))
} else{
	other_sample = setdiff(names(freq_two_cnts[1:2]), samp)
	pass_artif_other = pass_artif_snv[!(pass_artif_snv$Sample %in% c(samp, other_sample))
					& pass_artif_snv$pass == F, ]
	pass_artif_other$chr = gsub('chr', '', pass_artif_other$chr)
	recurrent_artif1 = names(which(table(pass_artif_other$id) >= 1))

	pass_artif_other = pass_artif_snv[pass_artif_snv$Sample != samp & pass_artif_snv$pass == F, ]
	pass_artif_other$chr = gsub('chr', '', pass_artif_other$chr)
	recurrent_artif2 = names(which(table(pass_artif_other$id) >= 2))

	recurrent_artif = union(recurrent_artif1, recurrent_artif2)
}

filt_recur_snv = samp_passes_snv1[samp_passes_snv1$id %in% recurrent_artif,]

# Assign passes
filters_snv = do.call(rbind, list(filt_bulkgatk_snv, filt_bulkmp_snv, filt_cls_snv, filt_recur_snv))
samp_passes_snv2 = samp_passes_snv1[!(samp_passes_snv1$id %in% filters_snv$id),]
gatk_snv$pass = F
gatk_snv[id %in% samp_passes_snv2$id, pass := T]

# Recompute the mutation burden
scan2_obj_newburden = recompute_burden(obj = res,
				calls_pass_snv = samp_passes_snv2,
				calls_pass_indel = samp_passes_indel2,
				cigar_change = T, lysis_change = F,
				cov_change = F,
				cigar_cutoff_fc = 1.5,
				quantile1 = 0.3, quantile2 = 0.8)

# Get final passes
samp_passes_snv_final = scan2_obj_newburden@gatk[scan2_obj_newburden@gatk$pass == T & scan2_obj_newburden@gatk$muttype == 'snv', ]

samp_passes_indel_final = scan2_obj_newburden@gatk[scan2_obj_newburden@gatk$pass == T & scan2_obj_newburden@gatk$muttype == 'indel', ]

# Get calls filtered due to final threshold (coverage, preamplification (aka lysis), cigar filter)
filt_final_snv = samp_passes_snv2[!(samp_passes_snv2$id %in% samp_passes_snv_final$id),]
filt_final_indel = samp_passes_indel2[!(samp_passes_indel2$id %in% samp_passes_indel_final$id),]

# Combine filtered snvs
if(nrow(filt_bulkgatk_snv) > 0){filt_bulkgatk_snv$reason = 'in_bulk_gatk'}else{filt_bulkgatk_snv = NULL}
if(nrow(filt_bulkmp_snv) > 0){filt_bulkmp_snv$reason = 'in_bulk_mpileup'}else{filt_bulkmp_snv = NULL}
if(nrow(filt_recur_snv) > 0){filt_recur_snv$reason = 'recurrent_snv'}else{filt_recur_snv = NULL}
if(nrow(filt_cls_snv) > 0){filt_cls_snv$reason = 'nearby_other_snvs'}else{filt_cls_snv = NULL}
if(nrow(filt_final_snv) > 0){filt_final_snv$reason = 'final_threshold'}else{filt_final_snv = NULL}

# Combine filtered indels
if(nrow(filt_bulkgatk_indel) > 0){filt_bulkgatk_indel$reason = 'in_bulk_gatk'}else{filt_bulkgatk_indel = NULL}
if(nrow(filt_morebulk_indel) > 0){filt_morebulk_indel$reason = 'more_indels_in_bulk'}else{filt_morebulk_indel = NULL}
if(nrow(filt_recur_indel) > 0){filt_recur_indel$reason = 'recurrent_indel'}else{filt_recur_indel = NULL}
if(nrow(filt_cls_indel) > 0){filt_cls_indel$reason = 'nearby_other_indels'}else{filt_cls_indel = NULL}
if(nrow(filt_final_indel) > 0){filt_final_indel$reason = 'final_threshold'}else{filt_final_indel = NULL}

comb_filt_indel = do.call(rbind, list(filt_bulkgatk_indel, filt_morebulk_indel,
				filt_recur_indel, filt_cls_indel, filt_final_indel))

comb_filt_snv = do.call(rbind, list(filt_bulkgatk_snv, filt_bulkmp_snv, filt_recur_snv,
				filt_cls_snv, filt_final_snv))

if(!is.null(comb_filt_indel)){
	comb_filt_indel$Sample = samp
	comb_filt_indel$sc.gt = comb_filt_indel[, .SD, .SDcols = samp]
}

if(!is.null(comb_filt_snv)){
	comb_filt_snv$Sample = samp
	comb_filt_snv$sc.gt = comb_filt_snv[, .SD, .SDcols = samp]
}

if(!is.null(samp_passes_snv_final)){
	samp_passes_snv_final$Sample = samp
	samp_passes_snv_final$sc.gt = samp_passes_snv_final[, .SD, .SDcols = samp]
}

if(!is.null(samp_passes_indel_final)){
	samp_passes_indel_final$Sample = samp
	samp_passes_indel_final$sc.gt = samp_passes_indel_final[, .SD, .SDcols = samp]
}

burden_df_snv = scan2_obj_newburden@mutburden$snv %>%
		add_column(quantile = c('Q1', 'Q2', 'Q3')) %>% 
		add_column(Sample = samp) %>%
		add_column(muttype = 'snv')

burden_df_indel = scan2_obj_newburden@mutburden$indel %>%
		add_column(quantile = c('Q1', 'Q2', 'Q3')) %>% 
		add_column(Sample = samp) %>%
		add_column(muttype = 'indel')

burden_df = rbind(burden_df_snv, burden_df_indel)

# Save
resL = list(burden_df, samp_passes_snv_final, samp_passes_indel_final, comb_filt_snv, comb_filt_indel)
names(resL) = c('Burden', 'Passes_snv', 'Passes_indel', 'Filtered_snv', 'Filtered_indel')
saveRDS(resL, paste0(outdir, '/', samp, '_final_burdens_and_calls.rds'))
saveRDS(scan2_obj_newburden, paste0(outdir, '/', samp, '_scan2_object_refine.rds'))




