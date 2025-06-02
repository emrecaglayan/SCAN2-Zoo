#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)


if (!(length(args) %in% 2:3)){
    stop("usage: save_passes_and_artifacts.R config.yaml outdir [n.cores]")
}

yaml <- args[1]
outdir <- args[2]
n.cores <- 1
if (length(args) == 3){
    n.cores <- as.integer(args[3])
}

library(scan2)
library(parallel)
library(data.table)

# Read panel config and extract params
yaml = yaml::read_yaml(yaml)

# Find and confirm scan2 outputs
scandirs = gsub('gatk*', 'call_mutations', yaml$sc_gvcfs)
scan2_outsL = sapply(scandirs, function(x){list.files(x, recursive = T, pattern = 'scan2_object_update.rda', full.names = T)})
scan2_outs = as.character(unlist(scan2_outsL))

if(length(scan2_outsL) > length(scan2_outs)){
	print("Some of the final scan2 objects are missing from some cells in the panel:")
	print(paste0("Number of scan2 objects: ", length(scan2_outs)))
	print(paste0("Number of cells in the panel: ", length(scan2_outsL)))
}

# For now, source the scripts that extract info from scan2 objects
scripts = paste0(yaml$scripts, '/UTILITY_FUNCTIONS.R')
scripts = gsub('SCAN2_BOWEN_NP', 'SCAN2_BOWEN_NP_TMP1', scripts)
source(scripts)

# Combine passes
sc_samps = gsub('.*mutations/|/scan2.*', '', scan2_outs)
passes = scan2_to_pass(scan2_outs[1:8], sc_samps[1:8], ncores = n.cores)

# Combine artifacts (also includes passes)
pass_artif = scan2_to_lenient_pass(scan2_files = scan2_outs[1:8],
					samps = sc_samps[1:8],
					ncores = n.cores,
					lysis_cutoff = 1)

# Save passes and artifacts
setcolorder(passes, c(setdiff(names(passes), "Sample"), "Sample"))
setcolorder(pass_artif, c(setdiff(names(pass_artif), "Sample"), "Sample"))

saveRDS(passes, paste0(outdir, '/AllPasses.RDS'))
saveRDS(pass_artif, paste0(outdir, '/AllPassesArtifacts.RDS'))



