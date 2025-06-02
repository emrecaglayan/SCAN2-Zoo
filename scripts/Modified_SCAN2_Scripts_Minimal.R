
## TO REMOVE PHASING REQUIREMENT ##

make.integrated.table2 = function(mmq60.tab, mmq1.tab, bulk.sample, sc.samples, 
    genome, snv.min.bulk.dp, indel.min.bulk.dp, snv.max.bulk.alt = 0, 
    snv.max.bulk.af = 0, indel.max.bulk.alt = 0, indel.max.bulk.af = 0, 
    panel = NULL, grs = genome.string.to.tiling(genome, tilewidth = 1e+07, 
        group = "auto"), legacy = FALSE, quiet = TRUE, report.mem = FALSE) 
{
    cat("Starting integrated table pipeline on", length(grs), 
        "chunks.\n")
    cat("Parallelizing with", future::nbrOfWorkers(), "cores.\n")
	progressr::with_progress({
        p <- progressr::progressor(along = 1:length(grs))
        p(amount = 0, class = "sticky", perfcheck(print.header = TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i, ]
            pc <- perfcheck(paste("read and annotate raw data", 
                i), {
                gatk <- read.tabix.data(path = mmq60.tab, region = gr, 
                  quiet = quiet, colClasses = list(character = "chr"))
                sitewide <- gatk[, 1:6]
                samplespecific <- gatk[, -(1:6)]
                annotate.gatk.counts(gatk.meta = sitewide, gatk = samplespecific, 
                  bulk.sample = bulk.sample, sc.samples = sc.samples, 
                  legacy = legacy, quiet = quiet)
                  
               # Modification: make add.mutsig=F for species that are not supported by SigProfiler
               if(!( grs@seqinfo@genome[1] %in% c("hs37d5", "hg38", "mm10","rn6")) ){
			annotate.gatk(gatk = sitewide, genome.string = genome, add.mutsig = F)
		} else{
               annotate.gatk(gatk = sitewide, genome.string = genome, add.mutsig = T)
               }
               
                annotate.gatk.lowmq(sitewide, path = mmq1.tab, 
                  bulk = bulk.sample, region = gr, quiet = quiet)
                  
               # Modification: Skip annotating phasing information
                #annotate.gatk.phasing(sitewide, phasing.path = phased.vcf, 
                #  region = gr, quiet = quiet)
                
                annotate.gatk.panel(sitewide, panel.path = panel, 
                  region = gr, quiet = quiet)
                annotate.gatk.candidate.loci(sitewide, snv.min.bulk.dp = snv.min.bulk.dp, 
                  snv.max.bulk.alt = snv.max.bulk.alt, snv.max.bulk.af = snv.max.bulk.af, 
                  indel.min.bulk.dp = indel.min.bulk.dp, indel.max.bulk.alt = indel.max.bulk.alt, 
                  indel.max.bulk.af = indel.max.bulk.af, mode = ifelse(legacy, 
                    "legacy", "new"))
            }, report.mem = report.mem)
            p(class = "sticky", amount = 1, pc)
            cbind(sitewide, samplespecific)
        })
    })
    gatk <- rbindlist(xs)
    if (nrow(gatk[somatic.candidate == TRUE]) == 0) 
        stop("0 somatic candidates detected. SCAN2 requires somatic candidates to have 0 supporting reads in 
the matched bulk - perhaps your bulk is too closely related to your single cells?")

    # Modification: See the modified function for details
    resampling.details <- gatk.resample.phased.sites(gatk)
    
    # Order columns as: 1 to somatic candidate index, resampled.training.site, bulk genotype and sc genotypes
    somcandInd <- which(colnames(gatk) == 'somatic.candidate')
    bulkInd <- which(colnames(gatk) == bulk.sample)
    bulkInds = bulkInd:(bulkInd+2)
    rtsInd = which(colnames(gatk) == 'resampled.training.site')
    the_rest = setdiff(1:ncol(gatk), c(1:(somcandInd), rtsInd, bulkInds))
    data.table::setcolorder(gatk, neworder = c(1:(somcandInd), rtsInd, bulkInds, the_rest))

    list(gatk = gatk, resampling.details = resampling.details)
}




gatk.resample.phased.sites2 = function (gatk, M = 20, seed = 0) 
{

	# Modification: Instead of properly phased snps, hsnps are chosen based on their af and depth using binomial distribution

	# Remove homozygous
	tokeep1 = which(!(gatk$balt == 0 | gatk$bulk.gt == '1/1' | gatk$bulk.gt == '0/0'))
	#gatk_filt <- gatk[!(balt == 0 | bulk.gt == '1/1' | bulk.gt == '0/0')]

	require(parallel)
	tokeep1_perc <- quantile(tokeep1, probs = seq(0, 1, by = 0.025))
	indexes <- as.numeric( sapply(tokeep1_perc, function(p) which.min(abs(tokeep1 - p))) )
	indexes = tokeep1[indexes]
	
	resL = mclapply(tokeep1, mc.cores = M, function(x){
			if(x %in% indexes){print(paste0('Choosing likely hSNPs without phasing: ', round(100 * (x / max(tokeep1))), '%'))}
			res = binom.confint(gatk[x,bulk.dp] * 0.5, gatk[x,bulk.dp], 0.99, method = 'bayes')
			interval_check = gatk[x,bulk.af] > res$lower & gatk[x,bulk.af] < res$upper
			nan_check = is.nan(res$lower) | is.nan(res$upper)
			!(nan_check) & interval_check
			})
	tokeep2 = tokeep1[which(unlist(resL))]

	# Subset the data frame to include only high-quality SNPs
	gatk[, hetSNP := .I %in% tokeep2]

	ret <- list()
	for (mt in c("snv", "indel")) {
		aux.data <- resample.germline(sites = gatk[somatic.candidate == T & muttype == mt],
						hsnps = gatk[hetSNP == T & muttype == mt], M = M, seed = seed)
		
		#is_present = paste0(gatk$chr, gatk$pos) %in% paste0(hqSNP$chr, hqSNP$pos)
		gatk[hetSNP == T & muttype == mt, `:=`(resampled.training.site, aux.data$selection$keep)]
		ret[[mt]] <- aux.data
	}
	gatk[is.na(resampled.training.site), `:=`(resampled.training.site, FALSE)]
	gatk[, hetSNP := NULL]
	return(ret)
}



convert_to_na <- function(arg) {
  if (arg == "NA") {
    return(NA)
  } else {
    return(arg)
  }
}



abmodel.fit.one.chrom2 = function(path, sc.sample, chrom, genome.seqinfo, n.chunks = 1, 
    n.logp.samples.per.chunk = 20000, hsnp.tilesize = 100, n.tiles = 250, 
    refine.n.steps = 4, refine.top.n = 50, logp.samples.per.step = 20000, 
    alim = c(-7, 2), blim = c(2, 4), clim = c(-7, 2), dlim = c(2, 
        6)) 
{
    chrom <- as.character(chrom)
    if (!(chrom %in% seqnames(genome.seqinfo))) {
        stop(paste0("invalid chromosome name '", chrom, "'\n"))
    }
    region <- as(genome.seqinfo, "GRanges")[chrom, ]
    hsnps <- read.training.hsnps(path = path, sample.id = sc.sample, 
        region = region)
    hsnps <- abmodel.downsample.hsnps(hsnps, hsnp.tilesize = hsnp.tilesize, 
        n.tiles = n.tiles, verbose = TRUE)
        
    # Assign max(scalt, scref) to scalt
    maxes = pmax(hsnps$scalt, hsnps$scref)
    mins = pmin(hsnps$scalt, hsnps$scref)
    hsnps$scalt = maxes
    hsnps$scref = mins

    # Assign max(scalt, scref) to phased.hap1
    hsnps$phased.hap1 = maxes
    hsnps$phased.hap2 = mins


    refine.records <- list()
    for (i in 1:refine.n.steps) {
        cat("Chrom:", chrom, "refinement step", i, "n.chunks", 
            n.chunks, "\n")
        refine.record <- abmodel.refine.parameter.space(hsnps, 
            n.chunks = n.chunks, n.logp.samples.per.chunk = n.logp.samples.per.chunk, 
            top.n = refine.top.n, hsnp.tilesize = hsnp.tilesize, 
            alim = alim, blim = blim, clim = clim, dlim = dlim)
        refine.records[[i]] <- refine.record
        alim <- refine.record$new.param.space[, "a"]
        blim <- refine.record$new.param.space[, "b"]
        clim <- refine.record$new.param.space[, "c"]
        dlim <- refine.record$new.param.space[, "d"]
    }
    list(n.chunks = n.chunks, n.logp.samples.per.chunk = n.logp.samples.per.chunk, 
        n.tiles = n.tiles, hsnp.tilesize = hsnp.tilesize, n.steps = refine.n.steps, 
        top.n = refine.top.n, refine.records = refine.records, 
        fit = refine.records[[refine.n.steps]]$logp.samples[1, 
            , drop = FALSE])
}


read.and.annotate.integrated.table2 <- function(path, sample.id, region=NULL, quiet=FALSE, ncores = 10) {

	tr <- read.integrated.table.1sample(path, sample.id, region = region, quiet = quiet)
	
	if(nrow(tr) == 0){
		return(tr)
	}
	
	setindex(tr, resampled.training.site)
	tr[, `:=`(dp, scalt + scref)]
	tr[, `:=`(af, scalt/dp)]
	data.table::setkey(tr, chr, pos, refnt, altnt)
	setindex(tr, muttype)
	sc.gt <- tr[[sample.id]]

	# Modification: Instead of properly phased snps, hsnps are chosen based on their af and depth using binomial distribution

	# Remove homozygous
	tokeep1 = which(!(tr$balt == 0 | tr$bulk.gt == '1/1' | tr$bulk.gt == '0/0'))

	require(parallel)
	tokeep1_perc <- quantile(tokeep1, probs = seq(0, 1, by = 0.025))
	indexes <- as.numeric( sapply(tokeep1_perc, function(p) which.min(abs(tokeep1 - p))) )
	indexes = tokeep1[indexes]
	
	resL = mclapply(tokeep1, mc.cores = ncores, function(x){
			if(x %in% indexes){print(paste0('Choosing likely hSNPs without phasing: ', round(100 * (x / max(tokeep1))), '%'))}
			res = binom.confint(tr[x,bulk.dp] * 0.5, tr[x,bulk.dp], 0.99, method = 'bayes')
			interval_check = tr[x,bulk.af] > res$lower & tr[x,bulk.af] < res$upper
			nan_check = is.nan(res$lower) | is.nan(res$upper)
			!(nan_check) & interval_check
			})
	tokeep2 = tokeep1[which(unlist(resL))]

	# Subset the data frame to include only high-quality SNPs
	tr[, hetSNP := .I %in% tokeep2]

	tr[, `:=`(training.site, hetSNP == T & sc.gt != "./." & bulk.gt != "./.")]
	tr[, hetSNP := NULL]
	# Set phased.hap1 and phased.hap2 in modified abmodel.fit.one.chrom function as max(scalt,scref)
	#tr[, `:=`(c("phased.hap1", "phased.hap2"), list(ifelse(phased.gt == "0|1", scref, scalt), ifelse(phased.gt == "0|1", scalt, scref)))]
	tr

}

read.integrated.table.1sample2 = function(path, sample.id, region = NULL, quiet = FALSE) 
{
    read.table.1sample(path = path, sample.id = sample.id, region = region, 
        meta.cols = integrated.table.meta.cols, quiet = quiet)
}

read.table.1sample2 <- function(path, sample.id, meta.cols, region=NULL, quiet=FALSE){
    if (!quiet) cat("Importing GATK table..\n")

	if(quiet == F){
		print(meta.cols)
	}
	
    n.meta.cols <- length(meta.cols)

    # Step 1: just get the header and detect the columns corresponding to sample.id
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- read.tabix.header(tf)
    close(tf)
    col.strings <- strsplit(header, '\t')[[1]]
    tot.cols <- length(col.strings)
    sample.idx <- which(col.strings == sample.id)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            if (i <= n.meta.cols) {
                cat(sprintf("    (%d)", i), col.strings[i], '\n')
            } else if (any(i == sample.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[single cell]\n')
            }
        }
    }
    
    # Step 2: really read the tables in, but only the relevant columns
    #cols.to.read <- rep("NULL", tot.cols)
    #cols.to.read[1:n.meta.cols] <- meta.cols

    # Read 3 columns for the single cell, 3 columns for bulk
    #cols.to.read[sample.idx + 0:2] <- c('character', 'integer', 'integer')

    # remove colClasses option
    gatk <- read.tabix.data(path=path, region=region, header=header, quiet=quiet)

    new.sample.idx <- which(colnames(gatk) == sample.id)
    colnames(gatk)[new.sample.idx+1:2] <- c('scref', 'scalt')

    # Rearrange columns so that the single cell triplet is first, then bulk triplet
    cols.to.keep <- col.strings[1:n.meta.cols]
    cols.to.keep <- c(cols.to.keep, sample.id, c('scref', 'scalt'))

    gatk <- gatk[,..cols.to.keep]
    gatk$chr = as.character(gatk$chr)

    gatk
}

read.integrated.table = function(object, path, quiet = FALSE) 
{
    object@gatk <- read.and.annotate.integrated.table(path = path, 
        sample.id = object@single.cell, region = object@region, 
        quiet = quiet)
    object@integrated.table.path <- path

    id <- object@gatk[, paste(chr, pos, refnt, altnt)]
    dupid <- duplicated(id)
    if (sum(dupid)) {
        stop(paste("found", sum(dupid), "duplicate (chr,pos,refnt,altnt) instances in integrated table. This is likely a bug, please report it."))
    }
    object
}


compute.ab.given.sites.and.training.data2 = function(sites, training.hsnps, ab.fits, quiet = FALSE) 
{
    chroms <- unique(sites$chr)

    if(nrow(training.hsnps) == 0){
    	return(NULL)
    }

    # Assign max(scalt, scref) to scalt
    maxes = pmax(training.hsnps$scalt, training.hsnps$scref)
    mins = pmin(training.hsnps$scalt, training.hsnps$scref)
    training.hsnps$scalt = maxes
    training.hsnps$scref = mins

    # Assign max(scalt, scref) to phased.hap1
    training.hsnps$phased.hap1 = maxes
    training.hsnps$phased.hap2 = mins

    do.work <- function(chrom, sites, ab.fit, hsnps) {
        if (!quiet) 
            cat(sprintf("inferring AB for %d sites on %s:%d-%d\n", 
                nrow(sites), chrom, as.integer(min(sites$pos)), 
                as.integer(max(sites$pos))))
              
        time.elapsed <- system.time(z <- infer.gp1(ssnvs = sites, 
            fit = ab.fit, hsnps = hsnps, flank = 1e+05, verbose = !quiet))
        if (!quiet) 
            print(time.elapsed)
        z
    }
    if (length(chroms) == 1) {
                print(head(training.hsnps))
        ab <- do.work(chrom = chroms, sites = sites, ab.fit = ab.fits[chroms, 
            , drop = FALSE], hsnps = training.hsnps)
    }
    else {
        ab <- do.call(rbind, lapply(chroms, function(chrom) {
            hsnps <- training.hsnps[chr == chrom]
            ab.fit <- ab.fits[chrom, , drop = FALSE]
            do.work(chrom = chrom, sites = sites, ab.fit = ab.fit, 
                hsnps = hsnps)
        }))
    }
    if (is.null(ab)) {
        data.frame(gp.mu = numeric(0), gp.sd = numeric(0))
    }
    else {
        ab
    }
}



run.pipeline2 = function(sc.sample, bulk.sample, int.tab, abfits, sccigars, 
    bulkcigars, trainingcigars, dptab, genome, grs = genome.string.to.tiling(genome, 
        tilewidth = 1e+07, group = "auto"), target.fdr = 0.01, 
    config.yaml = NULL, legacy = FALSE, report.mem = TRUE, verbose = TRUE) 
{
	cat("Starting chunked SCAN2 pipeline on", length(grs), "chunks\n")
	cat("Setting OpenBLAS corecount to 1. This prevents multithreaded matrix multiplication in chunks wher
	e it is undesired.\n")
	RhpcBLASctl::blas_set_num_threads(1)
	cat("Parallelizing with", future::nbrOfWorkers(), "cores\n")
	cat("Detailed chunk schedule:\n")
	cat(sprintf("%7s %5s %10s %10s\n", "Chunk", "Chr", "Start", 
	"End"))


	for (i in 1:length(grs)){
	cat(sprintf("%7d %5s %10d %10d\n", i, as.character(seqnames(grs)[i]), 
	    start(grs)[i], end(grs)[i]))
	}


	# Select regions that have some coverage for hetSNPs in the single cell WGS
	indsL = list()
	for(i in 1:length(grs)){
		gr <- grs[i, ]
		tr <- read.integrated.table.1sample(path = int.tab, sample.id = sc.sample, region = gr, quiet = T)
		tokeep1 = which(!(tr$balt == 0 | tr$bulk.gt == '1/1' | tr$bulk.gt == '0/0'))
		tr = tr[tokeep1]

		if(nrow(tr) > 0){
			indsL[[i]] = i
		}
		if(i%%20 == 0){print(paste0('Getting genomic regions to call mutations: ', i))}
	}


	grs_cov = unlist(indsL)

	progressr::with_progress({
	p <- progressr::progressor(along = grs_cov)
	p(amount = 0, class = "sticky", perfcheck(print.header = TRUE))
	xs <- future.apply::future_lapply(grs_cov, function(i) {
	    gr <- grs[i, ]
	    pc <- perfcheck(paste("make.scan", i), gt <- make.scan(single.cell = sc.sample, 
		bulk = bulk.sample, genome = genome, region = gr), 
		report.mem = report.mem)
	    p(class = "sticky", amount = 0, pc)
	    if (!is.null(config.yaml)) {
		gt <- add.static.filter.params(gt, config.path = config.yaml)
	    }
	    else {
		warning("config.yaml not specified, using R function default static filter parameters, which may not match the scan2 tool")
		gt <- add.static.filter.params(gt, muttype = "snv")
		gt <- add.static.filter.params(gt, muttype = "indel")
	    }
	    pc <- perfcheck(paste("read.integrated.table", i), 
		gt <- read.integrated.table(gt, path = int.tab, 
		  quiet = !verbose), report.mem = report.mem)
		gt@gatk$chr = as.character(gt@gatk$chr)
	    p(class = "sticky", amount = 0, pc)

	    gt@gatk$chr <- as.character(gt@gatk$chr)

	    pc <- perfcheck(paste("add.ab.fits", i), gt <- add.ab.fits(gt, 
		path = abfits), report.mem = report.mem)
	    p(class = "sticky", amount = 0, pc)
	    pc <- perfcheck(paste("compute.ab.estimates", i), 
		gt <- compute.ab.estimates(gt, quiet = !verbose), 
		report.mem = report.mem)
	    p(class = "sticky", amount = 0, pc)
	    pc <- perfcheck(paste("add.cigar.data", i), gt <- add.cigar.data(gt, 
		sccigars, bulkcigars, quiet = !verbose), report.mem = report.mem)
	    p(class = "sticky", amount = 0, pc)
	    pc <- perfcheck(paste("compute.models", i), gt <- compute.models(gt, 
		verbose = verbose), report.mem = report.mem)
	    p(class = "sticky", amount = 0, pc)
	    pc <- perfcheck(paste("compute.excess.cigar.scores", 
		i), gt <- compute.excess.cigar.scores(object = gt, 
		path = trainingcigars, legacy = TRUE, quiet = !verbose), 
		report.mem = report.mem)
	    p(class = "sticky", amount = 0, pc)
	    pc <- perfcheck(paste("compute.static.filters", i), 
		gt <- compute.static.filters(gt, mode = ifelse(legacy, 
		  "legacy", "new")), report.mem = report.mem)
	    p(class = "sticky", amount = 0, pc)
	    p()
	    gt
	}, future.seed = 0)
	})
	cat("Chunked pipeline complete.\n")
	x <- do.call(concat, xs)
	x <- compute.fdr.prior.data(x, mode = ifelse(legacy, "legacy", 
	"new"), quiet = !verbose)
	x <- compute.fdr(x, mode = ifelse(legacy, "legacy", "new"), 
	quiet = !verbose)
	x <- call.mutations(x, target.fdr = target.fdr, quiet = !verbose)
	x <- add.depth.profile(x, depth.path = dptab)
	x <- compute.mutburden(x)
	x
}


## TO REMOVE dbSNP REQUIREMENT ##
compute.filter.reasons2 <- function(gatk, target.fdr = o@call.mutations$target.fdr) 
{
    m <- gatk[, .(pass = !pass, abc.test, lysis.test = !is.na(lysis.fdr) & 
        lysis.fdr <= target.fdr, mda.test = !is.na(mda.fdr) & 
        mda.fdr <= target.fdr, cigar.id.test, cigar.hs.test, 
        lowmq.test, dp.test, min.sc.alt.test, csf.test)]
    m[is.na(m)] <- FALSE
    gatk[, `:=`(filter.reasons, apply(!m, 1, function(row) paste(colnames(m)[row], 
        collapse = "&")))]
}

#parse.static.filter.params2 <- function(config) {
#    setNames(lapply(c('snv', 'indel'), function(mt) {
#        list(
#            min.sc.alt=config[[paste0(mt, '_min_sc_alt')]],
#            min.sc.dp=config[[paste0(mt, '_min_sc_dp')]],
#            max.bulk.alt=config[[paste0(mt, '_max_bulk_alt')]],
#            max.bulk.af=config[[paste0(mt, '_max_bulk_af')]],
#            max.bulk.binom.prob=config[[paste0(mt, '_max_bulk_binom_prob')]],
#            min.bulk.dp=config[[paste0(mt, '_min_bulk_dp')]],
#            cg.id.q=config[[paste0(mt, '_cigar_id_score_quantile_cutoff')]],
#            cg.hs.q=config[[paste0(mt, '_cigar_hs_score_quantile_cutoff')]],
#            # Some cross-sample panel parameters that are currently not configurable
#            # via the usual process. For now, use update.static.filter.params() to change.
#            panel.use=ifelse(mt == 'snv', FALSE, TRUE),    # apply the panel or not? 
#            # IMPORTANT!! allow.missing was historically FALSE.  It is only recommended
#            # for very large panels.
#            panel.allow.missing=TRUE,                   # if a site is missing, does it PASS?
#            panel.max.unique.donors=1,
#            panel.max.out=2
#        )
#    }), c('snv', 'indel'))
#}



setValidity("SCAN2", function(object) {
    if (length(object@genome.string) != 1)
        return("must provide exactly one genome name")

    if (!is.null(object@genome.seqinfo)) {
        if (class(object@genome.seqinfo) != 'Seqinfo')
            stop('@genome.seqinfo must be of class Seqinfo')
    }

    if (length(object@single.cell) != 1)
        return("must provide exactly one single cell sample name")

    if (length(object@bulk) != 1)
        return("must provide exactly one bulk sample name")

    if (!is.null(object@region)) {
        if (!is(object@region, 'GRanges')) 
            stop('"region" must be a GRanges object')
        # This isn't possible if we're to combine chunks for the final
        # pipeline steps. Need a check that allows either single region
        # GRanges or the expected full autosome set.
        #
        # Disabling for now, but this isn't safe in general.
        #if (length(object@region) != 1) {
           #stop('"region" must contain exactly 1 range')
        #}
    }

    if (!is.null(object@gatk)) {
        if (ncol(object@gatk) != 13) {
            return("@gatk must have 13 columns")
        }
        if (!all(colnames(object@gatk)[1:7] ==
                c('chr', 'pos', 'refnt', 'altnt', 'mq', 'mqrs'))) {
            return("@gatk is improperly formatted")
        }
    }

    if (!is.null(object@ab.fits)) {
        # XXX: should also check that relevant chromosomes are present
        if (!all(colnames(object@ab.fits == c('a', 'b', 'c', 'd', 'logp'))))
            return("@ab.fits is improperly formatted")
    }

    if (!is.null(object@ab.estimates)) {
        if (!all(c('ab', 'gp.mu', 'gp.sd') %in% colnames(object@gatk)))
            return("ab.estimates were not properly imported")
    }

    if (!is.null(object@mut.models)) {
        if (!all(c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta')
            %in% colnames(object@gatk)))
            return("mut.models were not properly calculated")
    }

    # XXX: add validator for rest of slots
    return(TRUE)
})


setGeneric("concat2", function(...) standardGeneric("concat2"))
setMethod("concat2", signature="SCAN2", function(...) {
	args <- list(...)
	if (length(args) == 0) stop("tried to concat() 0 objects")
	if (length(args) == 1) return(args[[1]])
	init <- args[[1]]
	ret <- make.scan(init@single.cell, init@bulk, genome = init@genome.string, region = reduce(do.call(c, lapply(args, function(a) a@region))))
	ret@gatk <- rbindlist(lapply(args, function(a) a@gatk))
	ensure.same <- function(l, slot.name, var.name, var.name2) {
	    if (missing(var.name) & !missing(var.name2)) 
		stop("var.name2 only valid if var.name also specified")
	    if (all(sapply(l, function(element) is.null(slot(element, slot.name))))) {
		return()
	    }
	    if (!missing(var.name)) {
		if (missing(var.name2)) {
		    if (!all(sapply(l, function(element) all(slot(l[[1]], slot.name)[[var.name]] == slot(element, slot.name)[[var.name]])))) 
		        stop(paste("list of SCAN2 chunks cannot be concatenated; slot", slot.name, "does not have consistent values for", var.name))
		}
		else {
		    if (!all(sapply(l, function(element) all(slot(l[[1]], slot.name)[[var.name]][[var.name2]] == slot(element, slot.name)[[var.name]][[var.name2]])))) 
		        stop(paste("list of SCAN2 chunks cannot be concatenated; slot", slot.name, "does not have consistent values for", var.name))
		}
	    }
	    else {
		if (!all(sapply(l, function(element) all(slot(l[[1]], slot.name) == slot(element, slot.name))))) 
		    stop(paste("list of SCAN2 chunks cannot be concatenated; slot", slot.name, "does not have consistent values"))
	    }
	}
	ensure.same(args, "genome.string")
	ensure.same(args, "single.cell")
	ensure.same(args, "bulk")
	for (mt in c("snv", "indel")) {
	    ensure.same(args, "static.filter.params", mt, "min.sc.alt")
	    ensure.same(args, "static.filter.params", mt, "min.sc.dp")
	    ensure.same(args, "static.filter.params", mt, "min.bulk.dp")
	    ensure.same(args, "static.filter.params", mt, "max.bulk.alt")
	    ensure.same(args, "static.filter.params", mt, "max.bulk.af")
	    ensure.same(args, "static.filter.params", mt, "cg.id.q")
	    ensure.same(args, "static.filter.params", mt, "cg.hs.q")
	    ensure.same(args, "static.filter.params", mt, "disable.hs.filter")
	}
	ret@static.filter.params <- init@static.filter.params
	ensure.same(args, "cigar.data", "sc.path")
	ensure.same(args, "cigar.data", "bulk.path")
	ret@cigar.data <- data.frame(sc.sites = sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$sc.sites))), bulk.sites = sum(sapply(args, function(a) ifelse(is.null(a@cigar.data), 0, a@cigar.data$bulk.sites))))
	for (mt in c("snv", "indel")) {
	    ensure.same(args, "excess.cigar.scores", mt, "legacy")
	    ensure.same(args, "excess.cigar.scores", mt, "null.sites")
	}
	ret@excess.cigar.scores <- init@excess.cigar.scores
	ret@excess.cigar.scores$snv$sites <- sum(sapply(args, function(a) ifelse(is.null(a@excess.cigar.scores$snv), 0, a@excess.cigar.scores$snv$sites)))
	ret@excess.cigar.scores$indel$sites <- sum(sapply(args, function(a) ifelse(is.null(a@excess.cigar.scores$indel), 0, a@excess.cigar.scores$indel$sites)))
	for (mt in c("snv", "indel")) {
	    ensure.same(args, "fdr.prior.data", mt, "bins")
	    ensure.same(args, "fdr.prior.data", mt, "max.dp")
	    ensure.same(args, "fdr.prior.data", mt, "candidates.used")
	    ensure.same(args, "fdr.prior.data", mt, "hsnps.used")
	    ensure.same(args, "fdr.prior.data", mt, "nt.tab")
	    ensure.same(args, "fdr.prior.data", mt, "na.tab")
	    ensure.same(args, "fdr.prior.data", mt, "mode")
	}
	ret@fdr.prior.data <- init@fdr.prior.data
	ret@fdr <- NULL
	if (any(sapply(args, function(a) !is.null(a@fdr)))) {
	    ret@fdr <- list()
	    for (mt in c("snv", "indel")) {
		ensure.same(args, "fdr", mt, "mode")
		ret@fdr[[mt]] <- data.frame(mode = init@fdr[[mt]]$mode, sites = sum(sapply(args, function(a) ifelse(is.null(a@fdr[[mt]]), 0, a@fdr[[mt]]$sites))))
	    }
	}
	if (any(sapply(args, function(a) any(as.matrix(init@ab.fits) != as.matrix(a@ab.fits))))) stop("@ab.fits must be identical for all concat() elements")
	ret@ab.fits <- init@ab.fits
	ret@ab.estimates <- data.frame(sites = sum(sapply(args, function(a) ifelse(is.null(a@ab.estimates), 0, a@ab.estimates$sites))))
	ret@mut.models <- data.frame(sites = sum(sapply(args, function(a) ifelse(is.null(a@mut.models), 0, a@mut.models$sites))))
	ret
})



setGeneric("compute.fdr.prior.data2", function(object, mode=c('legacy', 'new'), quiet=FALSE)
    standardGeneric("compute.fdr.prior.data2"))
setMethod("compute.fdr.prior.data2", "SCAN2", function(object, mode=c('legacy', 'new'), quiet=FALSE)
{
	check.slots(object, c("gatk", "static.filter.params"))
	check.chunked(object, "compute.fdr.prior must be called on a SCAN2 object containing all sites, not a chunked parallelized object")
	mode <- match.arg(mode)
	if (!("static.filter" %in% colnames(object@gatk))) stop("must compute static filters before running compute.fdr.prior.data")
	muttypes <- c("snv", "indel")
	object@fdr.prior.data <- c(setNames(lapply(muttypes, function(mt) {
	    if (mode == "legacy") {
		sfp <- object@static.filter.params[["snv"]]
		cand <- object@gatk[muttype == mt & balt <= sfp$max.bulk.alt & bulk.gt == "0/0" & scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp & (is.na(balt.lowmq) | balt.lowmq <= sfp$max.bulk.alt)]
		hets = object@gatk[muttype == mt & training.site == TRUE & scalt >= sfp$min.sc.alt]
	    } else {
		sfp <- object@static.filter.params[[mt]]
		cand <- object@gatk[muttype == mt & somatic.candidate == TRUE & scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp]
		hets <- object@gatk[muttype == mt & training.site == TRUE & scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp]
	    }
	    c(compute.fdr.prior.data.for.candidates(candidates = cand, hsnps = hets, random.seed = 0, quiet = quiet), mode = mode)
	}), muttypes), list(mode = mode))
	object
})


setGeneric("compute.fdr2", function(object, path, mode=c('legacy', 'new'), quiet=FALSE)
    standardGeneric("compute.fdr2"))
setMethod("compute.fdr2", "SCAN2", function(object, path, mode=c('legacy', 'new'), quiet=FALSE)
{
    mode <- match.arg(mode)
    check.slots(object, c("gatk", "ab.estimates", "mut.models"))
    if (!missing(path) & !is.null(slot(object, "fdr.prior.data"))) 
        stop("fdr.prior.data is already loaded; path to new fdr.prior.data would overwrite")
    if (!missing(path)) 
        object@fdr.prior.data <- get(load(path))
    check.slots(object, "fdr.prior.data")
    muttypes <- c("snv", "indel")
    object@fdr <- c(setNames(lapply(muttypes, function(mt) {
        object@gatk[muttype == mt, `:=`(c("nt", "na"), estimate.fdr.priors(.SD, 
            object@fdr.prior.data[[mt]], use.ghet.loo = FALSE))]
        object@gatk[muttype == mt & training.site == TRUE & scalt >= 
            object@static.filter.params[[mt]]$min.sc.alt, `:=`(c("nt", 
            "na"), estimate.fdr.priors(.SD, object@fdr.prior.data[[mt]], 
            use.ghet.loo = TRUE))]
        if (mode == "legacy") {
            sfp <- object@static.filter.params[["snv"]]
            object@gatk[muttype == mt & (resampled.training.site | 
                (balt == 0 & bulk.gt == "0/0" &
                  scalt >= sfp$min.sc.alt & dp >= sfp$min.sc.dp)), 
                `:=`(c("lysis.fdr", "mda.fdr"), compute.fdr.legacy(altreads = scalt, 
                  dp = dp, gp.mu = match.ab(af = af, gp.mu = gp.mu), 
                  gp.sd = gp.sd, nt = nt, na = na, verbose = !quiet))]
            sites <- nrow(object@gatk[muttype == mt & balt == 
                0 & bulk.gt == "0/0" & scalt >= 
                sfp$min.sc.alt & dp >= sfp$min.sc.dp])
        } else if (mode == "new") {
            object@gatk[muttype == mt, `:=`(c("lysis.fdr", "mda.fdr"), 
                list(lysis.pv * na/(lysis.pv * na + lysis.beta * 
                  nt), mda.pv * na/(mda.pv * na + mda.beta * 
                  nt)))]
            sites <- nrow(object@gatk[muttype == mt])
        }
        list(sites = sites)
    }), muttypes), list(mode = mode))
    object
})



#setGeneric("compute.static.filters2", function(object, mode=c('new', 'legacy'))
#        standardGeneric("compute.static.filters2"))
#setMethod("compute.static.filters2", "SCAN2", function(object, mode=c('new', 'legacy'))
#{
#    check.slots(object, c("gatk", "cigar.data", "mut.models"))
#    mode <- match.arg(mode)
#    object@static.filter.params$mode <- mode
#    for (mt in c("snv", "indel")) {
#        sfp <- object@static.filter.params[[mt]]
#        object@gatk[muttype == mt, `:=`(c("cigar.id.test", "cigar.hs.test", 
#            "lowmq.test", "dp.test", "abc.test", "min.sc.alt.test", 
#            "max.bulk.alt.test", "max.bulk.af.test", "csf.test"), list(id.score >= object@excess.cigar.scores[[mt]]$id.score.q, 
#            hs.score >= object@excess.cigar.scores[[mt]]$hs.score.q, 
#            is.na(balt.lowmq) | balt.lowmq <= sfp$max.bulk.alt, 
#            dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp, 
#            abc.pv > 0.05, scalt >= sfp$min.sc.alt, balt <= sfp$max.bulk.alt, 
#            bulk.af <= sfp$max.bulk.af, muttype == "snv" | unique.donors <= 
#                1 | max.out <= 2))]
#        if (mode == "legacy") {
#            object@gatk[muttype == mt, `:=`(c("cigar.id.test", 
#                "cigar.hs.test"), list(id.score > object@excess.cigar.scores[[mt]]$id.score.q, 
#                hs.score > object@excess.cigar.scores[[mt]]$hs.score.q))]
#        }
#        object@gatk[muttype == mt, `:=`(static.filter, cigar.id.test & 
#            cigar.hs.test & lowmq.test & dp.test & abc.test & 
#            min.sc.alt.test & max.bulk.alt.test & max.bulk.af.test & csf.test)]
#    }
#    object
#
#})




compute.static.filters99 = function(object, mode=c('new', 'legacy'))
{
    check.slots(object, c("gatk", "cigar.data", "mut.models"))
    mode <- match.arg(mode)
    object@static.filter.params$mode <- mode
    for (mt in c("snv", "indel")) {
        sfp <- object@static.filter.params[[mt]]
        object@gatk[muttype == mt, `:=`(c("cigar.id.test", "cigar.hs.test", 
            "lowmq.test", "dp.test", "abc.test", "min.sc.alt.test", 
            "max.bulk.alt.test", "max.bulk.af.test", "csf.test"), list(id.score >= object@excess.cigar.scores[[mt]]$id.score.q, 
            hs.score >= object@excess.cigar.scores[[mt]]$hs.score.q, 
            is.na(balt.lowmq) | balt.lowmq <= sfp$max.bulk.alt, 
            dp >= sfp$min.sc.dp & bulk.dp >= sfp$min.bulk.dp, 
            abc.pv > 0.05, scalt >= sfp$min.sc.alt, balt <= sfp$max.bulk.alt, 
            bulk.af <= sfp$max.bulk.af, muttype == "snv" | unique.donors <= 
                1 | max.out <= 2))]
        if (mode == "legacy") {
            object@gatk[muttype == mt, `:=`(c("cigar.id.test", 
                "cigar.hs.test"), list(id.score > object@excess.cigar.scores[[mt]]$id.score.q, 
                hs.score > object@excess.cigar.scores[[mt]]$hs.score.q))]
        }
        object@gatk[muttype == mt, `:=`(static.filter, cigar.id.test & 
            cigar.hs.test & lowmq.test & dp.test & abc.test & 
            min.sc.alt.test & max.bulk.alt.test & max.bulk.af.test & csf.test)]
    }
    object

}



#helper.build.info.string2 <- function(gatk) {
#    info <- sprintf("MSC=%s;BALT_LOWMQ=%s",
#        ifelse(is.na(gatk$mutsig), '.', gatk$mutsig),
#        ifelse(is.na(gatk$balt.lowmq), '.', as.character(gatk$balt.lowmq)))
#    info <- paste0(info, ifelse(gatk$training.site, ';TS', ''))
#    info
#}



#filter.gatk.and.nearby.hets2 <- function(object, flank=2500, quiet=FALSE) {
#    if (!quiet) cat("Retaining candidate sites passing static filters and nearby germline hets...\n")
#    # new: try to retain mosaic sites. it is important to continue to retain all sites
#    # that could be rescued by mutation signature.  the bulk.af filter
#    som.idxs <- object@gatk[, which(((somatic.candidate & abc.test) | (bulk.binom.prob <= 1e-6 & bulk.af < 0.25)) & dp.test & min.sc.alt.test & (muttype == 'snv' | csf.test))]
#    flt <- object@gatk[som.idxs, .(chr, pos)]
#    gflt <- flank(GRanges(seqnames=flt$chr, ranges=IRanges(start=flt$pos, width=1)),
#        both=TRUE, width=flank)
#
#    gr <- GRanges(seqnames=object@gatk$chr, ranges=IRanges(start=object@gatk$pos, width=1))
#    ols <- countOverlaps(gr, gflt)
#    fl.idxs <- object@gatk[, which(ols > 0 & training.site == TRUE)]
#
#    # KEEP:
#    #   * all somatics sites passing above filters
#    #   * all training sites within `flank` bp of an above somatic site
#    #   * all resampled training sites, regardless of position
#    idxs.to.keep <- sort(union(union(som.idxs, fl.idxs),
#              which(object@gatk$resampled.training.site)))
#    ret <- object@gatk[idxs.to.keep]
#    compress.dt(ret)
#}


# Wrappers and expected formats for read.table.1sample
gatk.meta.cols2 <- c(
    chr='character',
    pos='integer',
    refnt='character',
    altnt='character'
)

integrated.table.meta.cols2 <- c(
    chr='character',
    pos='integer',
    refnt='character',
    altnt='character',
    mq = 'numeric',
    mqrs = 'numeric',
    bulk.gt='character',
    bref='integer',
    balt='integer',
    bulk.dp='integer',
    bulk.af='numeric',
#    bulk.binom.prob='numeric',
    tref='integer',
    talt='integer',
#    tcells='integer',
    muttype='character',
    mutsig='character',
    balt.lowmq='integer',
    nalleles='integer',
    unique.donors='integer',
    unique.cells='integer',
    unique.bulks='integer',
    max.out='integer',
    sum.out='integer',
    sum.bulk='integer',
    somatic.candidate='logical',
#    training.site='logical',
    resampled.training.site='logical'
)


annotate.gatk.candidate.loci2 <- function(gatk, snv.min.bulk.dp, snv.max.bulk.alt, snv.max.bulk.af, indel.min.bulk.dp, indel.max.bulk.alt, indel.max.bulk.af, mode = c("new", "legacy"))
{
    mode <- match.arg(mode)
    if (mode == "new") {
        snv.allow.bulk.gt <- snv.max.bulk.alt > 0 | snv.max.bulk.af > 
            0
        indel.allow.bulk.gt <- indel.max.bulk.alt > 0 | indel.max.bulk.af > 
            0
        gatk[, `:=`(somatic.candidate, ((muttype == "snv" & balt <= 
            snv.max.bulk.alt & bulk.dp >= snv.min.bulk.dp & !is.na(bulk.af) & 
            bulk.af <= snv.max.bulk.af & (is.na(balt.lowmq) | 
            balt.lowmq <= snv.max.bulk.alt) & (snv.allow.bulk.gt | 
            bulk.gt == "0/0")) | (muttype == "indel" & balt <= 
            indel.max.bulk.alt & bulk.dp >= indel.min.bulk.dp & 
            !is.na(bulk.af) & bulk.af <= indel.max.bulk.af & 
            (is.na(balt.lowmq) | balt.lowmq <= indel.max.bulk.alt) & 
            (indel.allow.bulk.gt | bulk.gt == "0/0"))) & talt >= 2)]
    }
    else if (mode == "legacy") {
        gatk[, `:=`(somatic.candidate, ((muttype == "snv" & balt <= 
            snv.max.bulk.alt) | (muttype == "indel" & balt <= 
            indel.max.bulk.alt)) & bulk.gt == "0/0" & talt >= 2)]
    }
    else stop("unrecognized mode")
    gatk
}



join.phased.hsnps2 <- function(bulk.called.vcf, hsnps.vcf, genome, grs = genome.string.to.tiling(genome, 
    tilewidth = 1e+07, group = "auto"), quiet = TRUE, report.mem = FALSE)
{
    cat("Starting phased hSNP joining on", length(grs), "chunks.\n")
    cat("Parallelizing with", future::nbrOfWorkers(), "cores.\n")
    progressr::with_progress({
        p <- progressr::progressor(along = 1:length(grs))
        p(amount = 0, class = "sticky", perfcheck(print.header = TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i, ]
            pc <- perfcheck(paste("read vcfs and join", i), {
                gatk.dt <- read.tabix.data(path = bulk.called.vcf, 
                  region = gr, quiet = quiet)[, 1:8]
                colnames(gatk.dt)[1:8] <- c("chr", "pos",
                  "refnt", "altnt", "qual", "filter", "info", 
                  "format")
                if (nrow(gatk.dt) == 0) {
                  gatk.dt[, `:=`(phasedgt, NA)]
                }
                else {
                  hsnps.dt <- read.tabix.data(path = hsnps.vcf, 
                    region = gr, quiet = quiet, header = "chr\tpos\trefnt\taltnt\tqual\tfilter\tinfo\tformat", 
                    colClasses = list(character = c("refnt", 
                      "altnt", "qual", "filter", "info")))
                  hsnps.dt[, `:=`(qual, ".")]
                  hsnps.dt[, `:=`(filter, ".")]
                  hsnps.dt[, `:=`(info, ".")]
                  hsnps.dt[, `:=`(format, "GT")]
                  gatk.dt[, `:=`(id, paste(chr, pos, refnt, 
                    altnt))]
                  hsnps.dt[, `:=`(id, paste(chr, pos, 
                    refnt, altnt))]
                  gatk.dt <- merge(hsnps.dt, gatk.dt[, .(id)], 
                    by = "id", all = FALSE, sort = FALSE)
                  gatk.dt$id <- NULL
                }
            }, report.mem = report.mem)
            p(class = "sticky", amount = 1, pc)
            gatk.dt
        })
    })
    gatk <- rbindlist(xs)
    gatk
}


## AVOIDS CONNECTION ERROR ##

retry_genome_string_to_tiling <- function(genome, tilewidth = 1e+07, group = "auto", max_attempts = 5) {
  attempts <- 0
  success <- FALSE
  result <- NULL
  
  while (!success && attempts < max_attempts) {
    attempts <- attempts + 1
    tryCatch({
      result <- genome.string.to.tiling(genome, tilewidth = tilewidth, group = group)
      success <- TRUE
    }, error = function(e) {
      cat("Attempt", attempts, "failed: ", e$message, "\n")
    })
  }
  
  if (!success) {
    stop("Failed to run genome.string.to.tiling after ", max_attempts, " attempts")
  }
  
  return(result)
}


digest.depth.profile2 = function(path, sc.sample, bulk.sample, genome, clamp.dp = 500, 
    grs = NULL, quiet = TRUE, report.mem = TRUE)
{

    grs = retry_genome_string_to_tiling(genome)
    print(grs)
    cat("Digesting depth profile using", length(grs), "chunks.\n")
    cat("Parallelizing with", future::nbrOfWorkers(), "cores.\n")
    progressr::with_progress({
        p <- progressr::progressor(along = 1:length(grs))
        p(amount = 0, class = "sticky", perfcheck(print.header = TRUE))
        xs <- future.apply::future_lapply(1:length(grs), function(i) {
            gr <- grs[i, ]
            pc <- perfcheck(paste("digest.depth.2sample", i), 
                dptab <- digest.depth.2sample(path = path, sc.sample = sc.sample, 
                  bulk.sample = bulk.sample, clamp.dp = clamp.dp, 
                  region = gr, quiet = quiet), report.mem = report.mem)
            p(class = "sticky", amount = 1, pc)
            dptab
        })
    }, enable = TRUE)
    dptab <- Reduce(`+`, xs)
    list(dptab = dptab, clamp.dp = clamp.dp)
}



## TO SUUPORT NON-HUMAN PRIMATE SPECIES ##


genome.string.to.seqinfo.object2 = function(genome = c("hs37d5", "hg38", "mm10", "rheMac10", "panTro6", "calJac4")) 
{
    genome <- match.arg(genome)
    if (genome == "hs37d5") {
        return(GenomeInfoDb::Seqinfo(genome = "GRCh37.p13"))
    }
    else if (genome == "hg38") {
        return(GenomeInfoDb::Seqinfo(genome = "hg38"))
    }
    else if (genome == "mm10") {
        return(GenomeInfoDb::Seqinfo(genome = "mm10"))
    }
    else if (genome == "rheMac10") {
        return(GenomeInfoDb::Seqinfo(genome = "rheMac10"))
    }
    else if (genome == "panTro6") {
        return(GenomeInfoDb::Seqinfo(genome = "panTro6"))
    }
    else if (genome == "calJac4") {
        return(GenomeInfoDb::Seqinfo(genome = "calJac4"))
    }
    else {
        stop("unsupported genome string")
    }
}



genome.string.to.tiling2 = function(genome = c("hs37d5", "hg38", "mm10", "rheMac10", "panTro6", "calJac4"), tilewidth = 1e+07, 
    group = c("auto", "sex", "circular", "all")) 
{
    genome <- match.arg(genome)
    group <- match.arg(group)
    if (genome == "hs37d5") {
        species <- "Homo_sapiens"
	sqi <- genome.string.to.seqinfo.object(genome)
	chroms.to.tile <- GenomeInfoDb::extractSeqlevelsByGroup(species = species,
	style = seqlevelsStyle(sqi)[1], group = group)
    }
    else if (genome == "hg38") {
        species <- "Homo_sapiens"
	sqi <- genome.string.to.seqinfo.object(genome)
	chroms.to.tile <- GenomeInfoDb::extractSeqlevelsByGroup(species = species,
	style = seqlevelsStyle(sqi)[1], group = group)
    }
    else if (genome == "mm10") {
        species <- "Mus_musculus"
	sqi <- genome.string.to.seqinfo.object(genome)
	chroms.to.tile <- GenomeInfoDb::extractSeqlevelsByGroup(species = species,
	style = seqlevelsStyle(sqi)[1], group = group)
    }
    else if (genome == "rheMac10") {
        species <- "Macaca_mulatta"
        
	# Get all chromosomes and haplotypes, keep only primary chromosome assemblies
	sqi <- genome.string.to.seqinfo.object(genome)
	sqi_levels = seqlevels(sqi)
	sqi_levels_uniq = unique(gsub('_.*', '', sqi_levels))
	lastchar = substr(sqi_levels_uniq, nchar(sqi_levels_uniq), nchar(sqi_levels_uniq))
	chroms.to.tile = sqi_levels_uniq[grep("[0-9]", lastchar)]
    }
    else if (genome == "panTro6") {
        species <- "Pan_troglodytes"
        
	# Get all chromosomes and haplotypes, keep only primary chromosome assemblies
	sqi <- genome.string.to.seqinfo.object(genome)
	sqi_levels = seqlevels(sqi)
	sqi_levels_uniq = unique(gsub('_.*', '', sqi_levels))
	lastchar = substr(sqi_levels_uniq, nchar(sqi_levels_uniq), nchar(sqi_levels_uniq))
	chroms.to.tile = sqi_levels_uniq[grep("[0-9]", lastchar)]
	chroms.to.tile = c(chroms.to.tile[1], 'chr2A', 'chr2B', chroms.to.tile[2:length(chroms.to.tile)])
    }
    else if (genome == "calJac4") {
        species <- "Callithrix_jacchus"
        
	# Get all chromosomes and haplotypes, keep only primary chromosome assemblies
	sqi <- genome.string.to.seqinfo.object(genome)
	sqi_levels = seqlevels(sqi)
	sqi_levels_uniq = unique(gsub('_.*', '', sqi_levels))
	lastchar = substr(sqi_levels_uniq, nchar(sqi_levels_uniq), nchar(sqi_levels_uniq))
	chroms.to.tile = sqi_levels_uniq[grep("[0-9]", lastchar)]
    }
    else {
        stop("unsupported genoome string")
    }
    
        
    grs <- GenomicRanges::tileGenome(seqlengths = sqi[chroms.to.tile], 
        tilewidth = tilewidth, cut.last.tile.in.chrom = TRUE)
    grs
}


genome.string.to.bsgenome.object2 = function(genome = c("hs37d5", "hg38", "mm10", 'rheMac10', 'panTro6', 'calJac4'))                                 
{
    genome <- match.arg(genome)
    if (genome == "hs37d5") {
        require(BSgenome.Hsapiens.1000genomes.hs37d5)                          
        genome <- BSgenome.Hsapiens.1000genomes.hs37d5                         
    }
    else if (genome == "hg38") {
        require(BSgenome.Hsapiens.UCSC.hg38)                                   
        genome <- BSgenome.Hsapiens.UCSC.hg38                                  
    }
    else if (genome == "mm10") {
        require(BSgenome.Mmusculus.UCSC.mm10)                                  
        genome <- BSgenome.Mmusculus.UCSC.mm10                                 
    }
    else if (genome == "rheMac10") {
        require(BSgenome.Mmulatta.UCSC.rheMac10)                                  
        genome <- BSgenome.Mmulatta.UCSC.rheMac10                               
    }
    else if (genome == "panTro6") {
        require(BSgenome.Ptroglodytes.UCSC.panTro6)                                  
        genome <- BSgenome.Ptroglodytes.UCSC.panTro6                                 
    }
    else if (genome == "calJac4") {
        require(BSgenome.Cjacchus.UCSC.calJac4)                                  
        genome <- BSgenome.Cjacchus.UCSC.calJac4                                 
    }
    else {
        stop("unsupported genome string")                                      
    }
    genome
}

make.scan2 = function(single.cell, bulk, genome = c("hs37d5", "hg38", "mm10"), region = NULL) 
{
    genome <- match.arg(genome, choices = c("hs37d5", "hg38", "mm10", 'rheMac10', 'panTro6', 'calJac4'))
    new("SCAN2", single.cell = single.cell, bulk = bulk, genome.string = genome, 
        genome.seqinfo = genome.string.to.seqinfo.object(genome), 
        region = region, gatk = NULL, ab.fits = NULL, ab.estimates = NULL, 
        mut.models = NULL, cigar.data = NULL, excess.cigar.scores = NULL, 
        static.filter.params = NULL, fdr.prior.data = NULL, fdr = NULL, 
        call.mutations = NULL, mutburden = NULL, mutsig.rescue = NULL, 
        spatial.sensitivity = NULL)
}

                                                                                  
get.gbp.by.genome2 = function (object)
{
    if (object@genome.string == "hs37d5") {                                                             
        total <- 3137161264                                                                             
        chrx <- 155270560
        chry <- 59373566                                                                                
        chrm <- 16571                                                                                   
        return((total - chrx - chry - chrm) * 2/1e+09)                                                  
    }
    else if (object@genome.string == "hg38") {                                                          
        total <- 3209286105                                                                             
        chrx <- 156040895
        chry <- 57227415                                                                                
        chrm <- 16569                                                                                   
        return((total - chrx - chry - chrm) * 2/1e+09)                                                  
    }                                                                                                   
    else if (object@genome.string == "mm10") {                                                          
        total <- 2730871774                                                                             
        chrx <- 171031299
        chry <- 91744698
        chrm <- 16299                                                                                   
        return((total - chrx - chry - chrm) * 2/1e+09)                                                  
    }                                                                                            
    else if (object@genome.string == "rheMac10") {                                                          
        total <- 2853992905                                                                             
        chrx <- 153388924
        chry <- 11753682                                                                                
        chrm <- 16564                                                                                   
        return((total - chrx - chry - chrm) * 2/1e+09)
    }                                                                                             
    else if (object@genome.string == "calJac4") {                                                          
        total <- 2830217926                                                                             
        chrx <- 148168104
        chry <- 8228174                                                                                
        chrm <- 16499                                                                                   
        return((total - chrx - chry - chrm) * 2/1e+09)
    }                                                                                              
    else if (object@genome.string == "panTro6") {                                                          
        total <- 2832631814                                                                             
        chrx <- 151576176
        chry <- 26350515                                                                                
        chrm <- 16554                                                                                   
        return((total - chrx - chry - chrm) * 2/1e+09)
    }                                                                        
    else {                                                                                              
        warn(paste("gbp not yet implemented for genome", object@genome.string))                         
        warn("the mutation burden for this analysis is a placeholder!")                                 
        warn("DO NOT USE!")
        return(-1)
    }
}


mutsig.rescue2 = function (object.paths, add.muts, rescue.target.fdr = 0.01, artifact.sigs = list(snv = data(snv.artifact.signature.v3), indel = data(indel.artifact.signature.v1)), true.sig = NULL,  quiet = FALSE, report.mem = TRUE) {
    if (is.null(names(object.paths))) 
        stop("output RDAs must be specified in the `names()` of `object.paths`")
    already.exists <- sapply(names(object.paths), file.exists)
    if (any(already.exists)) 
        stop(paste("these files specified in `names()` of `object.paths` already exist, please delete them and rerun this pipeline:", 
            paste(names(object.paths)[already.exists], collapse = " ")))
    use.add.muts <- FALSE
    if (!missing(add.muts) & !is.null(add.muts)) {
        if (!("data.table" %in% class(add.muts)) | !("muttype" %in% 
            colnames(add.muts)) | !("mutsig" %in% colnames(add.muts))) 
            stop("add.muts must be a data.table with \"muttype\" and \"mutsig\" columns")
        use.add.muts <- TRUE
    }
    cat("Rescuing mutations by signature.\n")
    cat("Parallelizing with", future::nbrOfWorkers(), "cores.\n")
    cat("WARNING: this pipeline requires ~6 GB of RAM per core due to data.table behavior that does not allow assignment by reference without first copying the entire table.\n")
    cat("Step 1. Getting high confidence mutations from", length(object.paths), 
        "SCAN2 objects.\n")
    progressr::with_progress({
        p <- progressr::progressor(along = 1:length(object.paths))
        p(amount = 0, class = "sticky", perfcheck(print.header = TRUE))
        gatks <- future.apply::future_lapply(1:length(object.paths), 
            function(i) {
                pc <- perfcheck(paste("prepare.object", i), {
                  x <- get(load(object.paths[i]))
                  if (is.compressed(x)) 
                    x <- decompress(x)
                  ret <- reduce.table(x@gatk, target.fdr = x@call.mutations$target.fdr)
                }, report.mem = report.mem)
                p(class = "sticky", amount = 1, pc)
                ret
            })
    }, enable = TRUE)
    cat("Step 2. Building true mutation spectra.\n")
    muttypes <- c("snv", "indel")
    true.sigs <- setNames(lapply(muttypes, function(mt) {
        if (!is.null(true.sig)) {
            return(true.sig[[mt]])
        }
        else {
            mutsigs <- do.call(c, lapply(gatks, function(gatk) get.high.quality.mutations(gatk, 
                muttype = mt)$mutsig))
            if (use.add.muts) {
                extra <- add.muts[muttype == mt]$mutsig
                cat(mt, ":", length(extra), "mutations taken from outside sources for true signature creation (add.muts)\n")
                mutsigs <- c(mutsigs, extra)
            }
            if (mt == "snv") 
                true.sig <- as.spectrum(sbs96(mutsigs))
            if (mt == "indel") 
                true.sig <- as.spectrum(id83(mutsigs))
            cat(mt, ": created true signature from", length(mutsigs), 
                "high confidence mutations.\n")
            return(true.sig)
        }
    }), muttypes)
    rm(gatks)
    gc()
    cat("Step3. Rescuing mutations and writing out new SCAN2 object RDA files.\n")
    progressr::with_progress({
        p <- progressr::progressor(along = 1:length(object.paths))
        p(amount = 0, class = "sticky", perfcheck(print.header = TRUE))
        results <- future.apply::future_lapply(1:length(object.paths), 
            function(i) {
                pc <- perfcheck(paste("mutsig.rescue.one", i), 
                  {
                    x <- get(load(object.paths[i]))
                    if (is.compressed(x)) 
                      x <- decompress(x)
                    x@mutsig.rescue <- NULL
                    x@gatk$rescue <- FALSE
                    for (mt in muttypes) {
                      x@mutsig.rescue[[mt]] <- mutsig.rescue.one(x, 
                        muttype = mt, artifact.sig = get(artifact.sigs[[mt]]), 
                        true.sig = true.sigs[[mt]], rescue.target.fdr = rescue.target.fdr)
                    }
                  }, report.mem = report.mem)
                p(class = "sticky", amount = 1, pc)
                results <- x
                save(results, file = names(object.paths)[i], 
                  compress = FALSE)
                results@gatk$sample <- results@single.cell
                results@gatk$bulk.sample <- results@bulk
                colnames(results@gatk)[colnames(results@gatk) == 
                  results@single.cell] <- "scgt"
                ret <- list(sample = results@single.cell, muts = results@gatk[pass == 
                  TRUE | rescue == TRUE], sig.homogeneity.test = setNames(lapply(results@mutsig.rescue, 
                  function(msr) msr$sig.homogeneity.test), names(results@mutsig.rescue)))
                ret
            })
    }, enable = TRUE)
    shts <- do.call(rbind, lapply(c("snv", "indel"), function(mt) {
        sht <- data.table(sample = sapply(results, function(r) r$sample), 
            muttype = mt, sig.homogeneity.test = sapply(results, 
                function(r) r$sig.homogeneity.test[[mt]]))
        sht$bonferroni <- p.adjust(p = sht$sig.homogeneity.test, 
            method = "bonferroni")
        sht$holm <- p.adjust(p = sht$sig.homogeneity.test, method = "holm")
        sht$bh.fdr <- p.adjust(p = sht$sig.homogeneity.test, 
            method = "fdr")
        sht
    }))
    list(muts = do.call(rbind, lapply(results, function(r) r$muts)), 
        sig.homogeneity.tests = shts)
}


genome.to.spmgr.format2 = c("GRCh37", "GRCh38", "mm10", "rheMac10", "panTro6", "calJac4")
names(genome.to.spmgr.format2) = c('hs37d', 'hg38', 'mm10', 'rheMac10', 'panTro6', 'calJac4')

####
## REPLACE THE OLD FUNCTIONS WITH NEW ONES IN THE NAMESPACE
####

unlockBinding("integrated.table.meta.cols", as.environment("package:scan2"))
assign("integrated.table.meta.cols", integrated.table.meta.cols2, envir = as.environment("package:scan2"))
if (exists("integrated.table.meta.cols", envir = .GlobalEnv)) {
  assign("integrated.table.meta.cols", integrated.table.meta.cols2, envir = .GlobalEnv)
}

assignInNamespace("make.integrated.table", make.integrated.table2, ns = "scan2")
if (exists("make.integrated.table", envir = .GlobalEnv)) {
  assign("make.integrated.table", make.integrated.table2, envir = .GlobalEnv)
}

assignInNamespace("abmodel.fit.one.chrom", abmodel.fit.one.chrom2, ns = "scan2")
if (exists("abmodel.fit.one.chrom", envir = .GlobalEnv)) {
  assign("abmodel.fit.one.chrom", abmodel.fit.one.chrom2, envir = .GlobalEnv)
}

assignInNamespace("read.and.annotate.integrated.table", read.and.annotate.integrated.table2, ns = "scan2")
if (exists("read.and.annotate.integrated.table", envir = .GlobalEnv)) {
  assign("read.and.annotate.integrated.table", read.and.annotate.integrated.table2, envir = .GlobalEnv)
}

assignInNamespace("read.table.1sample", read.table.1sample2, ns = "scan2")
if (exists("read.table.1sample", envir = .GlobalEnv)) {
  assign("read.table.1sample", read.table.1sample2, envir = .GlobalEnv)
}

assignInNamespace("gatk.resample.phased.sites", gatk.resample.phased.sites2, ns = "scan2")
if (exists("gatk.resample.phased.sites", envir = .GlobalEnv)) {
  assign("gatk.resample.phased.sites", gatk.resample.phased.sites2, envir = .GlobalEnv)
}

assignInNamespace("compute.ab.given.sites.and.training.data", compute.ab.given.sites.and.training.data2, ns = "scan2")
if (exists("compute.ab.given.sites.and.training.data", envir = .GlobalEnv)) {
  assign("compute.ab.given.sites.and.training.data", compute.ab.given.sites.and.training.data2, envir = .GlobalEnv)
}

assignInNamespace("run.pipeline", run.pipeline2, ns = "scan2")
if (exists("run.pipeline", envir = .GlobalEnv)) {
  assign("run.pipeline", run.pipeline2, envir = .GlobalEnv)
}

assignInNamespace("compute.filter.reasons", compute.filter.reasons2, ns = "scan2")
if (exists("compute.filter.reasons", envir = .GlobalEnv)) {
  assign("compute.filter.reasons", compute.filter.reasons2, envir = .GlobalEnv)
}


#parse.static.filter.params = parse.static.filter.params2

assignInNamespace("concat", concat2, ns = "scan2")
if (exists("concat", envir = .GlobalEnv)) {
  assign("concat", concat2, envir = .GlobalEnv)
}

assignInNamespace("concat", concat2, ns = "scan2")
if (exists("concat", envir = .GlobalEnv)) {
  assign("concat", concat2, envir = .GlobalEnv)
}

assignInNamespace("compute.fdr.prior.data", compute.fdr.prior.data2, ns = "scan2")
if (exists("compute.fdr.prior.data", envir = .GlobalEnv)) {
  assign("compute.fdr.prior.data", compute.fdr.prior.data2, envir = .GlobalEnv)
}

assignInNamespace("compute.fdr", compute.fdr2, ns = "scan2")
if (exists("compute.fdr", envir = .GlobalEnv)) {
  assign("compute.fdr", compute.fdr2, envir = .GlobalEnv)
}

assignInNamespace("compute.static.filters", compute.static.filters99, ns = "scan2")
if (exists("compute.static.filters", envir = .GlobalEnv)) {
  assign("compute.static.filters", compute.static.filters99, envir = .GlobalEnv)
}

#assignInNamespace("helper.build.info.string", helper.build.info.string2, ns = "scan2")
#if (exists("helper.build.info.string", envir = .GlobalEnv)) {
#  assign("helper.build.info.string", helper.build.info.string2, envir = .GlobalEnv)
#}

#helper.build.info.string = helper.build.info.string2

#assignInNamespace("write.vcf", write.vcf2, ns = "scan2")
#if (exists("write.vcf", envir = .GlobalEnv)) {
#  assign("write.vcf", write.vcf2, envir = .GlobalEnv)
#}

#write.vcf = write.vcf2

#assignInNamespace("filter.gatk.and.nearby.hets", filter.gatk.and.nearby.hets2, ns = "scan2")
#if (exists("filter.gatk.and.nearby.hets", envir = .GlobalEnv)) {
#  assign("filter.gatk.and.nearby.hets", filter.gatk.and.nearby.hets2, envir = .GlobalEnv)
#}

assignInNamespace("gatk.meta.cols", gatk.meta.cols2, ns = "scan2")
if (exists("gatk.meta.cols", envir = .GlobalEnv)) {
  assign("gatk.meta.cols", gatk.meta.cols2, envir = .GlobalEnv)
}

assignInNamespace("annotate.gatk.candidate.loci", annotate.gatk.candidate.loci2, ns = "scan2")
if (exists("annotate.gatk.candidate.loci", envir = .GlobalEnv)) {
  assign("annotate.gatk.candidate.loci", annotate.gatk.candidate.loci2, envir = .GlobalEnv)
}

assignInNamespace("join.phased.hsnps", join.phased.hsnps2, ns = "scan2")
if (exists("join.phased.hsnps", envir = .GlobalEnv)) {
  assign("join.phased.hsnps", join.phased.hsnps2, envir = .GlobalEnv)
}

assignInNamespace("read.integrated.table.1sample", read.integrated.table.1sample2, ns = "scan2")
if (exists("read.integrated.table.1sample", envir = .GlobalEnv)) {
  assign("read.integrated.table.1sample", read.integrated.table.1sample2, envir = .GlobalEnv)
}

assignInNamespace("digest.depth.profile", digest.depth.profile2, ns = "scan2")
if (exists("digest.depth.profile", envir = .GlobalEnv)) {
  assign("digest.depth.profile", digest.depth.profile2, envir = .GlobalEnv)
}


assignInNamespace("genome.string.to.seqinfo.object", genome.string.to.seqinfo.object2, ns = "scan2")
if (exists("genome.string.to.seqinfo.object", envir = .GlobalEnv)) {
  assign("genome.string.to.seqinfo.object", genome.string.to.seqinfo.object2, envir = .GlobalEnv)
}


assignInNamespace("genome.string.to.tiling", genome.string.to.tiling2, ns = "scan2")
if (exists("genome.string.to.tiling", envir = .GlobalEnv)) {
  assign("genome.string.to.tiling", genome.string.to.tiling2, envir = .GlobalEnv)
}

assignInNamespace("genome.string.to.bsgenome.object", genome.string.to.bsgenome.object2, ns = "scan2")
if (exists("genome.string.to.bsgenome.object", envir = .GlobalEnv)) {
  assign("genome.string.to.bsgenome.object", genome.string.to.bsgenome.object2, envir = .GlobalEnv)
}


unlockBinding("genome.to.spmgr.format", as.environment("package:scan2"))
assign("genome.to.spmgr.format", genome.to.spmgr.format2, envir = as.environment("package:scan2"))
if (exists("genome.to.spmgr.format", envir = .GlobalEnv)) {
  assign("genome.to.spmgr.format", genome.to.spmgr.format2, envir = .GlobalEnv)
}

assignInNamespace("make.scan", make.scan2, ns = "scan2")
if (exists("make.scan", envir = .GlobalEnv)) {
  assign("make.scan", make.scan2, envir = .GlobalEnv)
}


assignInNamespace("get.gbp.by.genome", get.gbp.by.genome2, ns = "scan2")
if (exists("get.gbp.by.genome", envir = .GlobalEnv)) {
  assign("get.gbp.by.genome", get.gbp.by.genome2, envir = .GlobalEnv)
}









