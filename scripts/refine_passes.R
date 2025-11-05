#!/usr/bin/env Rscript

# Arguments from bash
print(commandArgs(trailingOnly = TRUE))

args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) == 3)){
    stop("usage: refine_passes.R config.yaml scan2_object_update.rda merge_mmq60.tab.gz")
}

# Get from bash
yaml_fn = args[1]
scan2_fn = args[2]
panel_fn = args[3]

# Extract from config file
config = yaml::read_yaml(yaml_fn)
samp_bulk_bam = config$bulk_bam
gb = config$genome
bsgenome = config$bsgenome
fastafile = config$ref
id_frac = config$indel_fraction_delta
metafile = config$makepanel_metadata
stringency = config$stringency
cnt_alt = config$cnt_alt
q1 = config$q1
q2 = config$q2

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
meta = read.csv(metafile)

dir.create('objects', showWarnings=F)

tmpdir = 'tmpdir'
dir.create(tmpdir, showWarnings=F)

# Load gatk
res = get(load(scan2_fn))
sc = res@single.cell
gatk = res@gatk

# Add new metric to gatk and save back in the object
gatk$sc_id_frac = (gatk$ID.cigars/gatk$dp.cigars) - (gatk$ID.cigars.bulk/gatk$dp.cigars.bulk)
gatk$sc_id_frac_afNorm = gatk$sc_id_frac / gatk$af
gatk$id = paste(gatk$chr, gatk$pos, gatk$refnt, gatk$altnt, sep = '_')
res@gatk = gatk

gatk_indel = gatk[gatk$muttype == 'indel',]
gatk_snv = gatk[gatk$muttype == 'snv',]

rm(gatk)
gc()


####
## INDELS
####

## Remove if in bulk ##
samp = basename(dirname(scan2_fn))
outdir = paste0('objects/', samp)
dir.create(outdir, showWarnings = F)

print(paste0('Processing sample for indel ', samp))

## Add the ones that are in the bulk up to 2 reads ##
samp_passes_indel1 = gatk_indel[cigar.id.test & cigar.hs.test & dp.test &
	abc.test & min.sc.alt.test & csf.test & bulk.gt == '0/0' &
	mda.fdr < 0.01 & lysis.fdr < 0.01 & (balt.lowmq <= 2 | balt <= 2),]
print(paste0('Total indel passes ', sum(gatk_indel$pass)))
print(paste0('Total indel passes with relaxed balt.lowmq ', nrow(samp_passes_indel1)))


if(nrow(samp_passes_indel1) == 0){
	print('No indels found, skipping indel regenotyping...')
	print('Setting the indel variables to empty data table to avoid errors')
	samp_passes_indel2 = samp_passes_indel1
	filters_indel = samp_passes_indel1
	filt_bulkgatk_indel = samp_passes_indel1
	filt_cls_indel = samp_passes_indel1
	filt_morebulk_indel = samp_passes_indel1
	filt_bulkcigar = samp_passes_indel1
} else{

	## GATK - Regenotype new passes, remove if present in bulk ##
	sites = paste0(samp_passes_indel1$chr, ':', samp_passes_indel1$pos, '-', samp_passes_indel1$pos)

	# Write sites to bed file
	sites = paste0(samp_passes_indel1$chr, ':', samp_passes_indel1$pos)
	intervals_file = paste0(outdir, "/indel_sites.bed")
	intervals = sub(":", "\t", sub("-", "\t", sites))
	writeLines(intervals, intervals_file)

	regen_ids_bulk = regenotype_gatk(intervals_file = intervals_file, samp = samp, outdir = outdir,
				fastafile = fastafile, bamfile = samp_bulk_bam, pref = c('bulkbam', 'indel'))

	filt_bulkgatk_indel = samp_passes_indel1[(samp_passes_indel1$id %in% regen_ids_bulk),]

	## Remove calls in sites with many indels and clippings ##
	# This filter essentially removes indel calls that have many indels in both single-cell and bulk
	# before local realignment. So these are the indel calls that remain specific to single-cell after
	# the local realignment. While these may be true calls, the indels are usually not common indel types
	# as they are large deletions of repetitive regions. So I have decided to remove them with this filter.

	indel_cigarpass = gatk_indel[sc_id_frac_afNorm > id_frac, id]
	filt_morebulk_indel = samp_passes_indel1[!(samp_passes_indel1$id %in% indel_cigarpass),]

	# Remove if multiple indels are nearby #
	filt_cls_indel = find_clustered_snvs(obj = res, inpass = samp_passes_indel1)

	# Remove any sites with more than 2 indels in bulk, indicating two fragments of indels
	filt_bulkcigar = samp_passes_indel1[ID.cigars.bulk > 2,]

	filters_indel = do.call(rbind, list(filt_bulkgatk_indel, filt_morebulk_indel, filt_cls_indel, filt_bulkcigar))
	samp_passes_indel2 = samp_passes_indel1[!(samp_passes_indel1$id %in% filters_indel$id),]
	gatk_indel$pass = F
	gatk_indel[id %in% samp_passes_indel2$id, pass := T]
}


####
## SNVs
####

print(paste0('Processing sample for sSNV ', samp))

## Add the ones that are in the bulk up to 2 reads ##
samp_passes_snv1 = gatk_snv[cigar.id.test & cigar.hs.test & dp.test &
	abc.test & min.sc.alt.test & bulk.gt == '0/0' &
	mda.fdr < 0.01 & lysis.fdr < 0.01 & (balt.lowmq <= 2 | balt <= 2),]

if(gb == 'hg19'){
	samp_passes_snv1$chr = gsub('chr', '', samp_passes_snv1$chr)
}

print(paste0('Total snv passes ', sum(gatk_snv$pass)))
print(paste0('Total snv passes with relaxed balt.lowmq ', nrow(samp_passes_snv1)))

## GATK - Regenotype new passes, remove if present in bulk ##

# Write sites to bed file
sites = paste0(samp_passes_snv1$chr, ':', samp_passes_snv1$pos)
intervals_file = paste0(outdir, "/snv_sites.bed")
intervals = sub(":", "\t", sub("-", "\t", sites))
writeLines(intervals, intervals_file)


regen_ids_bulk = regenotype_gatk(intervals_file = intervals_file, samp = samp, outdir = outdir,
			fastafile = fastafile, bamfile = samp_bulk_bam, pref = c('bulkbam', 'snv'))

filt_bulkgatk_snv = samp_passes_snv1[(samp_passes_snv1$id %in% regen_ids_bulk),]


## MPILEUP - Regenotype new passes, remove if present in bulk ##
# This filter will remove a sSNV call if at least one alt-allele supporting read is present in bulk
# either before or after local realignment of the sSNV sites in bulk

bulk_pass = get_snv_mpileup_pass(samp_passes_snv1 = samp_passes_snv1, outdir = outdir, samp_bulk_bam = samp_bulk_bam, tmpdir = tmpdir)
filt_bulkmp_snv = samp_passes_snv1[!(samp_passes_snv1$id %in% bulk_pass),]

## Remove if multiple sSNVs are nearby ##
filt_cls_snv = find_clustered_snvs(obj = res, inpass = samp_passes_snv1)


## Remove if recurrent other samples in the panel ##
panel = fread(panel_fn)
colnames(panel)[1] = 'chr'
filt_recurrent = get_nxp_artifacts(passes = samp_passes_snv1, sc = sc, meta = metafile, panel = panel, stringency = stringency, cnt_alt = cnt_alt, gb = gb, bsgenome = bsgenome)
filt_recurrent$id0 = NULL


# Assign passes
filters_snv = do.call(rbind, list(filt_bulkgatk_snv, filt_bulkmp_snv, filt_cls_snv, filt_recurrent))
samp_passes_snv2 = samp_passes_snv1[!(samp_passes_snv1$id %in% filters_snv$id),]
gatk_snv$pass = F
gatk_snv[id %in% samp_passes_snv2$id, pass := T]


# Let's go with a modest increase in cigar cutoff as default. But this can be made more stringent
# if there is a clear trend of burden estimate going down with greater cigar cutoff stringency.
# If that is the case, it would mean that called mutations are removed more than hetSNPs with increasing
# cigar cutoffs, indicating that increasing cigar stringency likely removes more false positives.
scan2_obj_newburden = recompute_burden(obj = res,
				calls_pass_snv = samp_passes_snv2,
				calls_pass_indel = samp_passes_indel2,
				cigar_change = T, lysis_change = T,
				cov_change = F, lysis_thresh = 0.05,
				cigar_cutoff_fc = 2,
				quantile1 = q1, quantile2 = q2)

# If burden estimate failed, there is a good chance this is due to low data quality and low coverage.
# The loop below increases depth quantile until the burden is estimated.
while ((is.na(scan2_obj_newburden@mutburden$snv[2, 'ncalls']) || 
        is.na(scan2_obj_newburden@mutburden$indel[2, 'ncalls'])) && 
       (q2 < 1)) {
	
	inc = 0.01
	q1 = q1 + inc
	q2 = q2 + inc
	scan2_obj_newburden = recompute_burden(obj = res,
					calls_pass_snv = samp_passes_snv2,
					calls_pass_indel = samp_passes_indel2,
					cigar_change = T, lysis_change = T,
					cov_change = F, lysis_thresh = 0.05,
					cigar_cutoff_fc = 2,
					quantile1 = q1, quantile2 = q2)
	print(paste0('Lower quantile increased to ', q1))
}


# Get final passes
samp_passes_snv_final = scan2_obj_newburden@gatk[scan2_obj_newburden@gatk$pass == T & scan2_obj_newburden@gatk$muttype == 'snv', ]

# Plot mutation profile
if(gb == 'hg19'){
	plot_mutprof_from_vcf(mutdf = samp_passes_snv_final, genome_build = gb,
	ref_genome = bsgenome,
	prefix = paste0('MP_', sc), ymax1 = 0.15,
	oldcols = c('chr', 'pos', 'refnt', 'altnt'),
	add_chr = T)
} else{
	plot_mutprof_from_vcf(mutdf = samp_passes_snv_final, genome_build = gb,
	ref_genome = bsgenome,
	prefix = paste0('MP_', sc), ymax1 = 0.15,
	oldcols = c('chr', 'pos', 'refnt', 'altnt'),
	add_chr = F)
}

samp_passes_indel_final = scan2_obj_newburden@gatk[scan2_obj_newburden@gatk$pass == T & scan2_obj_newburden@gatk$muttype == 'indel', ]

# Get calls filtered due to final threshold (coverage, preamplification (aka lysis), cigar filter)
filt_final_snv = samp_passes_snv2[!(samp_passes_snv2$id %in% samp_passes_snv_final$id),]
filt_final_indel = samp_passes_indel2[!(samp_passes_indel2$id %in% samp_passes_indel_final$id),]

# Combine filtered snvs
if(nrow(filt_bulkgatk_snv) > 0){filt_bulkgatk_snv$reason = 'in_bulk_gatk'}else{filt_bulkgatk_snv = NULL}
if(nrow(filt_bulkmp_snv) > 0){filt_bulkmp_snv$reason = 'in_bulk_mpileup'}else{filt_bulkmp_snv = NULL}
if(nrow(filt_cls_snv) > 0){filt_cls_snv$reason = 'nearby_other_snvs'}else{filt_cls_snv = NULL}
if(nrow(filt_recurrent) > 0){filt_recurrent$reason = 'recurrent'}else{filt_recurrent = NULL}
if(nrow(filt_final_snv) > 0){filt_final_snv$reason = 'final_threshold'}else{filt_final_snv = NULL}

# Combine filtered indels
if(nrow(filt_bulkgatk_indel) > 0){filt_bulkgatk_indel$reason = 'in_bulk_gatk'}else{filt_bulkgatk_indel = NULL}
if(nrow(filt_morebulk_indel) > 0){filt_morebulk_indel$reason = 'more_indels_in_bulk'}else{filt_morebulk_indel = NULL}
if(nrow(filt_bulkcigar) > 0){filt_bulkcigar$reason = 'high_indels_in_bulk'}else{filt_bulkcigar = NULL}
if(nrow(filt_cls_indel) > 0){filt_cls_indel$reason = 'nearby_other_indels'}else{filt_cls_indel = NULL}
if(nrow(filt_final_indel) > 0){filt_final_indel$reason = 'final_threshold'}else{filt_final_indel = NULL}

comb_filt_indel = do.call(rbind, list(filt_bulkgatk_indel, filt_morebulk_indel, filt_cls_indel, filt_bulkcigar, filt_final_indel))

comb_filt_snv = do.call(rbind, list(filt_bulkgatk_snv, filt_bulkmp_snv, filt_cls_snv, filt_recurrent, filt_final_snv))

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
save(scan2_obj_newburden, file = paste0(outdir, '/', samp, '_scan2_object_refine.rda'))




