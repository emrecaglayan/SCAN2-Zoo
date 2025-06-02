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
            snakemake@input['mmq60'],
            snakemake@input['mmq1'],
            #snakemake@input['phasing'],
            # if the user did not provide a panel, this entry is NULL and c() drops
            # NULL values silently. replace with NA and handle it later.
            ifelse(is.null(snakemake@params[['cross_sample_panel']]), NA,
                snakemake@params[['cross_sample_panel']]),
            snakemake@params['config_yaml'],
            snakemake@params['bulk_sample'],
            snakemake@params['snv_min_bulk_dp'],
            snakemake@params['snv_max_bulk_alt'],
            snakemake@params['snv_max_bulk_af'],
            snakemake@params['indel_min_bulk_dp'],
            snakemake@params['indel_max_bulk_alt'],
            snakemake@params['indel_max_bulk_af'],
            snakemake@params['genome'],
            snakemake@params['legacy'],
            snakemake@output['tab'],
            snakemake@output['tabgz'],
            snakemake@output['details'],
            snakemake@threads
        ))
        ret
    }
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)

# Load packages
suppressMessages(library(scan2))
suppressMessages(library(yaml))
suppressMessages(library(future))
suppressMessages(library(binom))

# Load the modified script
source('/home/ch253794/workdir/PROGRAMS/SCAN2_BOWEN_NP/scripts/Modified_SCAN2_Scripts_Minimal.R')

# Old variables are commented out
mmq60 <- args[1]
mmq1 <- args[2]
panel <- args[3]
config.yaml <- args[4]
bulk.sample <- args[5]
snv.min.bulk.dp <- as.integer(args[6])
snv.max.bulk.alt <- as.integer(args[7])
snv.max.bulk.af <- as.numeric(args[8])
indel.min.bulk.dp <- as.integer(args[9])
indel.max.bulk.alt <- as.integer(args[10])
indel.max.bulk.af <- as.numeric(args[11])
genome <- args[12]
legacy <- as.logical(args[13])
out.tab <- args[14]
out.tab.gz <- args[15]
out.rda <- args[16]

n.cores <- 1
if (length(args) == 17)
    n.cores <- as.integer(args[17])

for (f in c(out.tab, out.tab.gz, paste0(out.tab.gz, '.tbi'), out.rda)) {
    if (file.exists(f))
        stop(paste('output file', f, 'already exists, please delete it first'))
}

plan(multicore, workers=n.cores)

if (is.na(panel))
    print('Panel is null')
    panel <- NULL

y <- yaml::read_yaml(config.yaml)
sc.samples <- names(y$sc_bams)

# Remove the dbsnp column from tables
mmq60_r = data.table::fread(mmq60)
if(any(colnames(mmq60_r) == 'dbsnp')){
	mmq60_r[,which(colnames(mmq60_r) == 'dbsnp')] = NULL
	fwrite(mmq60_r, sep = '\t', gsub('.gz', '', mmq60))
	Rsamtools::bgzip(gsub('.gz', '', mmq60), dest = mmq60, overwrite = T)
	Rsamtools::indexTabix(mmq60, format = "vcf", skip = 1)
}

mmq1_r = data.table::fread(mmq1)
if(any(colnames(mmq1_r) == 'dbsnp')){
	mmq1_r[,which(colnames(mmq1_r) == 'dbsnp')] = NULL
	fwrite(mmq1_r, sep = '\t', gsub('.gz', '', mmq1))
	Rsamtools::bgzip(gsub('.gz', '', mmq1), dest = mmq1, overwrite = T)
	Rsamtools::indexTabix(mmq1, format = "vcf", skip = 1)
}

# Currently the chunking used here isn't configurable by user
results <- make.integrated.table(mmq60.tab=mmq60, mmq1.tab=mmq1, panel=panel,
    bulk.sample=bulk.sample, sc.samples=sc.samples,
    snv.min.bulk.dp=snv.min.bulk.dp,
    indel.min.bulk.dp=indel.min.bulk.dp,
    snv.max.bulk.alt=snv.max.bulk.alt,
    snv.max.bulk.af=snv.max.bulk.af,
    indel.max.bulk.alt=indel.max.bulk.alt,
    indel.max.bulk.af=indel.max.bulk.af,
    genome=genome, legacy=legacy)

inttab <- results$gatk
write.integrated.table(inttab=inttab, out.tab=out.tab, out.tab.gz=out.tab.gz)

resampling.details <- results$resampling.details
save(resampling.details, file=out.rda)

if ('snakemake' %in% ls()) {
    sink()
}
