#!/usr/bin/env Rscript

update_integrate_table <- function(
    tab.path, panel.path, out_tab.path, genome.string,
    grs=genome.string.to.tiling(genome.string, tilewidth=10e6, group='auto'),
    quiet=TRUE, report.mem=FALSE){
    cat('Starting update rda pipeline on', length(grs), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i,]
            pc <- perfcheck(paste('update results data',i), {
		tab <- read.tabix.data(path=tab.path, region=gr, quiet=quiet,colClasses=list(character='chr'))		
		panel <- read.tabix.data(path=panel.path, region=gr, quiet=quiet, colClasses=list(character='chr'))
		tab[panel, on=.(chr,pos),
			c('nalleles', 'unique.donors', 'unique.cells', 'unique.bulks', 'max.out', 'sum.out', 'sum.bulk') :=
			list(i.nalleles, i.unique.donors, i.unique.cells, i.unique.bulks, i.max.out, i.sum.out, i.sum.bulk)]	

            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)
            tab
        })
    })
    tab <- rbindlist(xs)

    write.integrated.table(inttab=tab, out.tab=out_tab.path, out.tab.gz=paste0(out_tab.path,'.gz')) 
}


args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 4)
    stop('update_integrate_table.R in.tab panel.tab.gz out.rda [n.cores]')
in.rda <- args[1]
out.rda <-args[2]
panel.path  <- args[3]
genome <- args[4]

n.cores <- 1
if (length(args) == 5)
    n.cores <- as.integer(args[5])

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

library(scan2)
library(future)
library(progressr)
plan(multicore, workers=n.cores)

update_integrate_table(in.rda, panel.path, out.rda, genome)

