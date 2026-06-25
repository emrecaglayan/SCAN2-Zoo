#!/usr/bin/env Rscript

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) {
        ret <- unlist(c(
            snakemake@input['tab'],
            snakemake@params['sc_sample'],
            snakemake@params['bulk_sample'],
            snakemake@params['genome'],
            snakemake@params['config_yaml'],
            snakemake@output['rda'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 5) {
    stop("usage: digest_depth.R joint_depth_matrix.tab.gz sc.sample bulk.sample genome out.rda [n.cores]")
}

path <- args[1]
sc.sample <- args[2]
bulk.sample <- args[3]
genome <- args[4]
config.yaml <- args[5]
out.rda <- args[6]

n.cores <- 1
if (length(args) == 7)
    n.cores <- as.integer(args[7])

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

suppressMessages(library(scan2.zoo))
suppressMessages(library(future))
suppressMessages(library(progressr))
plan(multicore, workers=n.cores)

y <- yaml::read_yaml(config.yaml)

# Get all chromosomes that are in the config file. These should all be diploid chromosomes
scaffolds = unique(gsub(':.*', '', y$gatk_regions))
scaffolds = scaffolds[order(as.numeric(gsub("[^0-9]", "", scaffolds)))]
grs = genome.string.to.tiling(genome, tilewidth = 1e+07, scaffolds = scaffolds)

progressr::with_progress({
    # handler_newline causes alot of printing, but it's log-friendly
    progressr::handlers(progressr::handler_newline())
    results <- digest.depth.profile(path=path, sc.sample=sc.sample, bulk.sample=bulk.sample,
        grs=grs)
}, enable=TRUE)

dptab <- results$dptab
clamp.dp <- results$clamp.dp
save(dptab, clamp.dp, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
