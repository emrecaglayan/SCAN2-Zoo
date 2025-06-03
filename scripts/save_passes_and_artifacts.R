#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)


if (!(length(args) %in% 2:3)){
    stop("usage: save_passes_and_artifacts.R config.yaml [n.cores]")
}

yaml_fn <- args[1]
outdir <- 'artifacts_panel'
dir.create(outdir, showWarnings=F)

n.cores <- 1
if (length(args) == 3){
    n.cores <- as.integer(args[3])
}

print(n.cores)

suppressPackageStartupMessages(suppressWarnings(library(scan2)))
suppressPackageStartupMessages(suppressWarnings(library(parallel)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))

# Read panel config and extract params
config = yaml::read_yaml(yaml_fn)
metafile = config$makepanel_metadata

# Find and confirm scan2 outputs
scandirs = gsub('gatk*', 'call_mutations', config$sc_gvcfs)
scan2_outsL = sapply(scandirs, function(x){list.files(x, recursive = T, pattern = 'scan2_object_update.rda', full.names = T)})
scan2_outs = as.character(unlist(scan2_outsL))

if(length(scan2_outsL) > length(scan2_outs)){
	print("Some of the final scan2 objects are missing from some cells in the panel:")
	print(paste0("Number of scan2 objects: ", length(scan2_outs)))
	print(paste0("Number of cells in the panel: ", length(scan2_outsL)))
}

# For now, source the scripts that extract info from scan2 objects
scripts = paste0(config$scripts, '/UTILITY_FUNCTIONS.R')
scripts = gsub('SCAN2_BOWEN_NP', 'SCAN2_BOWEN_NP_TMP1', scripts)
source(scripts)

# Combine passes
sc_samps = gsub('.*mutations/|/scan2.*', '', scan2_outs)
passes = scan2_to_pass(scan2_outs, sc_samps, ncores = n.cores)

# Combine artifacts (also includes passes)
pass_artif = scan2_to_lenient_pass(scan2_files = scan2_outs,
					samps = sc_samps,
					ncores = n.cores,
					lysis_cutoff = 1)

# Combine with metadata
meta = read.csv(metafile)
meta$Sample = meta$sample
meta$sample = NULL

passes = merge(passes, meta, by = 'Sample')
pass_artif = merge(pass_artif, meta, by = 'Sample')


# Save passes and artifacts
setcolorder(passes, c(setdiff(names(passes), "Sample"), "Sample"))
setcolorder(pass_artif, c(setdiff(names(pass_artif), "Sample"), "Sample"))

saveRDS(passes, paste0(outdir, '/AllPasses.RDS'))
saveRDS(pass_artif, paste0(outdir, '/AllPassesArtifacts.RDS'))



