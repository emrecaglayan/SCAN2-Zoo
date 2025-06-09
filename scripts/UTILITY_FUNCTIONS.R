
compute_mutburden_new = function(object, gbp.per.genome = get.gbp.by.genome(object), 
    quiet = FALSE, q1 = 0.4, q2 = 0.9) 
{
    check.slots(object, c("call.mutations", "depth.profile"))
    muttypes <- c("snv", "indel")
    object@mutburden <- setNames(lapply(muttypes, function(mt) {

        pre.geno.burden <- object@fdr.prior.data[[mt]]$burden[2]
        sfp <- object@static.filter.params[[mt]]
        g <- object@gatk[resampled.training.site == TRUE & muttype == 
            mt]
        dptab <- object@depth.profile$dptab
        dptab <- dptab[1:min(max(g$dp) + 1, nrow(dptab)), ]
        if (nrow(g) < 100) {
            warning(paste("only", nrow(g), "resampled germline", 
                mt, "sites were detected; aborting genome-wide extrapolation. Typical whole-genome experiments include ~10-100,000 germline sites"))
            ret <- data.frame(ncalls = NA, callable.sens = NA, 
                callable.bp = NA)[c(1, 1, 1), ]
        }
        else {
            s <- object@gatk[pass == TRUE & muttype == mt]
            q = 4
            qstouse <- c(1, 2, 4)
            qbreaks <- quantile(g$dp, prob = c(0, q1, (q1+q2)/2, q2, 1.0))
            if (length(unique(qbreaks)) != q + 1) {
                ncalls <- rep(NA, 3)
                callable.sens <- rep(NA, 3)
                rowqs <- rep(NA, 3)
                callable.bp <- rep(NA, 3)
                warning(paste("could not derive unique breakpoints for quartiles of sequencing depth at germline hSNPs.  this usually indicates that sequencing depth is heavily skewed toward low depths (typically DP=0)\ngot qbreaks = ", 
                  deparse(qbreaks)))
            }
            else {
                s$dpq <- cut(s$dp, qbreaks, include.lowest = T, 
                  labels = F)
                s$dpq[s$dpq == 3] <- 2
                g$dpq <- cut(g$dp, qbreaks, include.lowest = T, 
                  labels = F)
                g$dpq[g$dpq == 3] <- 2
                rowqs <- cut(0:(nrow(dptab) - 1), qbreaks, include.lowest = T, 
                  labels = F)
                rowqs[rowqs == 3] <- 2
                s <- s[dpq %in% qstouse]
                g <- g[dpq %in% qstouse]
                ncalls <- sapply(qstouse, function(q) sum(s[dpq == 
                  q]$pass, na.rm = TRUE))
                callable.bp <- sapply(split(dptab[, -(1:sfp$min.bulk.dp)], 
                  rowqs), sum)
                callable.sens <- sapply(qstouse, function(q) mean(g[bulk.dp >= 
                  sfp$min.bulk.dp & dpq == q & resampled.training.site == 
                  TRUE]$training.pass, na.rm = TRUE))
            }
            ret <- data.frame(ncalls = ncalls, callable.sens = callable.sens, 
                callable.bp = callable.bp)
        }
        ret$callable.burden <- ret$ncalls/ret$callable.sens
        ret$rate.per.gb <- ret$callable.burden/ret$callable.bp * 
            1e+09/2
        ret$burden <- ret$rate.per.gb * gbp.per.genome
        ret$somatic.sens <- ret$ncalls/ret$burden
        ret$pre.genotyping.burden <- pre.geno.burden
        ret$unsupported.filters <- sfp$max.bulk.alt > 0 | sfp$max.bulk.af > 0
        if (any(ret$unsupported.filters)) {
            warning("mutation burdens must be extrapolated without clonal mutations (i.e., max.bulk.alt=0 and max.bulk.af=0)! burdens will be estimated, but they are invalid")
        }
        ret
    }), muttypes)
    object
}


recompute_burden = function(obj, calls_pass_snv, calls_pass_indel, cigar_change = F, lysis_change = F, lysis_thresh = 0.01, cigar_cutoff_fc = 1, cov_change = F, cov_thresh = 25, quantile1 = 0.25, quantile2 = 0.75){

	# Extract gatk
	raw_gatk = obj@gatk

	# Give unique id to sites
	raw_gatk[, id := paste0(chr, '_', pos, '_', refnt, '_', altnt)]
	calls_pass_snv$id = paste0(calls_pass_snv$chr, '_', calls_pass_snv$pos, '_', calls_pass_snv$refnt, '_', calls_pass_snv$altnt)
	calls_pass_indel$id = paste0(calls_pass_indel$chr, '_', calls_pass_indel$pos, '_', calls_pass_indel$refnt, '_', calls_pass_indel$altnt)

	# Assign new mutations
	raw_gatk$pass = F
	raw_gatk[(id %in% calls_pass_snv$id) | (id %in% calls_pass_indel$id), pass := T]

	if(cigar_change == T){

	for(mtype in c('snv', 'indel')){

		if(mtype == 'snv'){
			# Extract score cutoffs
			id_score_q = obj@excess.cigar.scores$snv$id.score.q
			hs_score_q = obj@excess.cigar.scores$snv$hs.score.q

			# Assign more strict cigar passes
			raw_gatk = raw_gatk[, pass_strictCigar := muttype == 'snv' &
					id.score > id_score_q*cigar_cutoff_fc &
					hs.score > hs_score_q*cigar_cutoff_fc, ]
			
			# Remove somatic mutation calls and training.passes if they fail the new cigar filter
			raw_gatk[muttype == 'snv' & pass == T & pass_strictCigar == F, pass := F]
			raw_gatk[muttype == 'snv' & training.pass == T & pass_strictCigar == F, training.pass := F]

		} else if(mtype == 'indel'){
			# Extract score cutoffs
			id_score_q = obj@excess.cigar.scores$indel$id.score.q
			hs_score_q = obj@excess.cigar.scores$indel$hs.score.q

			# Assign more strict cigar passes
			raw_gatk = raw_gatk[, pass_strictCigar := muttype == 'indel' &
					id.score > id_score_q*cigar_cutoff_fc &
					hs.score > hs_score_q*cigar_cutoff_fc, ]
			
			# Remove somatic mutation calls and training.passes if they fail the new cigar filter
			raw_gatk[muttype == 'indel' & pass == T & pass_strictCigar == F, pass := F]
			raw_gatk[muttype == 'indel' & training.pass == T & pass_strictCigar == F, training.pass := F]
		}
	}
	}

	if(lysis_change == T){
		raw_gatk[pass == T & lysis.pv >= lysis_thresh, pass := F]
		raw_gatk[training.pass == T & lysis.pv >= lysis_thresh, training.pass := F]
	}

	if(cov_change == T){
		raw_gatk[pass == T & dp <= cov_thresh, pass := F]
		raw_gatk[training.pass == T & dp <= cov_thresh, training.pass := F]
	}

	obj@gatk = raw_gatk
	old_mb = obj@mutburden

	obj = compute_mutburden_new(object = obj,
				gbp.per.genome = get.gbp.by.genome(obj),
				quiet = F, q1 = quantile1, q2 = quantile2)

	return(obj)
}



plot_mutprof_from_mutmat = function(mut_mat, genome_build, ref_genome,
				prefix = '', hg = 7, wd1 = 14, wd2 = 10,
				add_chr = F, ymax1 = 0.2, ymax2 = 1, plot = T,
				oldcols = c('chr', 'pos', 'refnt', 'altnt')){

	setnames(mutdf, old = oldcols, new = c("chr","start","REF","ALT"), skip_absent=TRUE)

	if(add_chr == T){mutdf$chr = paste0('chr', mutdf$chr)}

	mutdf$end = mutdf$start
	mutgr = makeGRangesFromDataFrame(mutdf, keep.extra.columns=T)
	GenomeInfoDb::genome(mutgr) = genome_build
	
	mut_mat <- mut_matrix(vcf_list = mutgr, ref_genome = ref_genome)

	if(plot == T){

		pdf(paste(prefix, 'MutationSignatures_TripleNucleotide.pdf', sep = '_'), width = wd1, height = hg)
		print( plot_96_profile2(mut_mat, plot_size_scale = 2, condensed = T, ymax = ymax1) )
		dev.off()
	}

	type_occurrences <- mut_type_occurrences(mutgr, ref_genome)
	
	# Combine two specific groups (e.g., "G1" and "G2")
	df = as.data.frame(t(type_occurrences))
	df$ctgroup = rownames(df)
	df <- df %>%
	  mutate(group_combined = ifelse(ctgroup %in% c("C>T at CpG", "C>T other"), "C>T", ctgroup))
	df = df[df$ctgroup != 'C>T',]
	df$relcount = df$My_sample / sum(df$My_sample)

	# Plot
	cols = c('skyblue', 'black', 'pink', 'red', 'gray', 'darkolivegreen3', 'mistyrose')

	if(plot == T){

		pdf(paste(prefix, 'MutationSignatures_SingleNucleotide.pdf', sep = '_'), width = wd2, height = hg)
		print(	

		ggplot(df, aes(x = group_combined, y = relcount, fill = ctgroup)) +
		geom_bar(stat = "identity", position = "stack") +
		theme_classic() +
		xlab('') + ylab('Relative Contribution') +
		scale_fill_manual(values = cols) +
		theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
		ylim(0, ymax2) +
		rotate_x_text(90)

		)
		dev.off()
	}

	return(mut_mat)
}




ppvenn = function(set1, set2, setnames){

	suppressPackageStartupMessages(require(VennDiagram))
	size1 = length(set1)
	size2 = length(set2)
	size3 = length(intersect(set1, set2))
	fn = paste0(paste(setnames, collapse = '_'), '_Overlap.pdf')

	pdf(fn, width = 8)
	print( draw.pairwise.venn(area1 = size2,
		area2 = size1,
		cross.area = size3, print.mode = c('raw', 'percent'),
		fill = c('orange', 'darkgreen'), col = c('orange', 'darkgreen'),
		category = c(setnames[1], setnames[2]), fontface = 'bold', cex = 2,
		cat.cex = 2, cat.fontface = 'bold', margin = 0.05) )
	dev.off()

}


find_clustered_snvs = function(obj = res, extsize = 5000, inpass){

	artpass = obj@gatk[obj@gatk$static.filter & obj@gatk$lysis.fdr < 0.9,]
	indL = list()
	for(k in 1:nrow(inpass)){

		snv = inpass[k, ]
		st = snv$pos - extsize
		end = snv$pos + extsize
		snvs = artpass[chr == snv$chr & pos > st & pos < end,]
		if(nrow(snvs) > nrow(snv)){indL[[k]] = k}
	}
	return(inpass[unlist(indL),])
}

remove_bulk_mpileup = function(passes_file, bulk_bam_file, outfile, altnt = 'altnt', refnt = 'refnt', pos = 'pos', chr = 'chr'){

	tmpL = lapply(1:nrow(passes_file), function(x){

		if(passes_file[x, 'muttype'] == 'snv'){

		# Extract the snv and its basecalls
		site = passes_file[x,]
		alt = passes_file[x, altnt]
		ref_basecall = extract_basecall_and_read(bam = bulk_bam_file, outfile = outfile,
								chr = as.character(site[,chr]),
								pos = as.character(site[,pos]))
		}
		if(x%%50==0){print(x)}

		ref_basecall$V3 = gsub("(.)\\1+", "\\1", ref_basecall$V3)
		sum(ref_basecall$V3 %in% alt) / nrow(ref_basecall)
	})

	# Fail if it returns null
	tmpL = lapply(tmpL, function(x){ifelse(length(x)==0, 0, x)})
	passes_file_bulkFilt = passes_file[!(unlist(tmpL) > 0),]
	return(passes_file_bulkFilt)
}

scan2_to_lenient_pass = function(scan2_files, samps, ncores = 2, lysis_cutoff = 0.5){

	suppressPackageStartupMessages(require(dplyr))
	suppressPackageStartupMessages(require(parallel))
	# Load lysis fails from this database
	gatkLysisL = mclapply(1:length(scan2_files), mc.cores = ncores, function(x){
			obj = get(load(scan2_files[x]))
			gatk = obj@gatk
			gatk$Sample = samps[x]
			gc()
			print(paste0('Loaded cell ', x))
			gatk[static.filter == TRUE & lysis.fdr < lysis_cutoff & abc.test == TRUE]
	})

	cmnCols = Reduce(intersect, lapply(gatkLysisL, function(x){colnames(x)}))
	gatkLysis = lapply(gatkLysisL, function(x){x[, ..cmnCols]}) %>% do.call(rbind, .)

	gatkLysis$id = paste0(gatkLysis$chr, '_', gatkLysis$pos, '_', gatkLysis$refnt, '_', gatkLysis$altnt)

	return(gatkLysis)
}

add_mutsigs = function(mutdf, bsgenome, extension = 1, gb){

	suppressPackageStartupMessages(require(MutationalPatterns))
	suppressPackageStartupMessages(require(GenomicRanges))

	if(!grepl('chr', mutdf$chr[1])){
		mutdf$chr = paste0('chr', mutdf$chr)
	}
	setnames(mutdf, old = c('chr', 'pos', 'refnt', 'altnt'), new = c("chr","start","REF","ALT"), skip_absent=TRUE)
	mutdf$end = mutdf$start
	mutdf$pos = NULL
	mutdf2 = makeGRangesFromDataFrame(mutdf, keep.extra.columns = T)
	GenomeInfoDb::genome(mutdf2) = gb
	type_context = type_context(mutdf2, bsgenome, extension = 1)
	mutsigs = paste0(type_context[[2]], ':', type_context[[1]])
	mutdf[mutdf$muttype == 'snv', 'mutsig'] = mutsigs
	mutdf$end = NULL
	setnames(mutdf, new = c('chr', 'pos', 'refnt', 'altnt'), old = c("chr","start","REF","ALT"), skip_absent=TRUE)
	return(mutdf)
}

safe_read_table <- function(file, ...) {
  tryCatch(
    {
      read.table(file, ...)
    },
    error = function(e) {
      message("Error reading file: ", file, "\nReturning NULL.")
      return(NULL)
    }
  )
}

mutrate_comp_plot = function(mmu_df1, hsa_df1, sample_1, sample_2, annot_pos_x = 0.15, annot_pos_y = 0.7, palette = c('blue', 'red')){

	mmu_df1$species = sample_1
	hsa_df1$species = sample_2
	toplot = rbind(mmu_df1, hsa_df1)

	p = ggscatter(toplot, x = 'age', y = 'burden', size = 2, color = 'species', add = 'reg.line',
	add.params = list(fill = "species"), conf.int = T, palette = palette) +
	stat_cor(aes(color = species), method = "pearson", label.x.npc = "left", label.y.npc = "top") +
	xlab('Age') + ylab('Mutation Burden') +
	theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20),
		axis.title = element_text(size=20)) +
		rotate_x_text(90)


	# Mutation rates comparison
	library(lme4)
	lm_mmu = lmer(burden ~ age + (1|subject), mmu_df1)
	lm_hsa = lmer(burden ~ age + (1|subject), hsa_df1)

	mmu_age = summary(lm_mmu)$coefficients['age', 'Estimate'] %>% round(., digits = 2)
	hsa_age = summary(lm_hsa)$coefficients['age', 'Estimate'] %>% round(., digits = 2)

	mmu_std = summary(lm_mmu)$coefficients['age', 'Std. Error'] %>% round(., digits = 2)
	hsa_std = summary(lm_hsa)$coefficients['age', 'Std. Error'] %>% round(., digits = 2)


	library(cowplot)
	custom_text <- ggdraw() +
	  draw_plot(p) +
	  draw_text(paste0("Mutation rates\n", sample_1, ": ", mmu_age, ", +- ", mmu_std,
	  			"\n", sample_2, ": ", hsa_age, ", +- ", hsa_std),  fontface = 'bold',
		    x = annot_pos_x, y = annot_pos_y, hjust = 0, size = 14)


	pdf(paste0(sample_1, '_', sample_2, '_NeuN_MutationRate.pdf'), width = 10)
	print(custom_text)
	dev.off()
}


plot_scan2_ab = function(scan2_fls, samps, ncores = 4, prefix, rand_n = 500000){

	suppressPackageStartupMessages(require(parallel))
	print('Getting hetSNPs from SCAN2 object')

	tmpL = mclapply(1:length(scan2_fls), mc.cores = ncores, function(x){
		res = get(load(scan2_fls[x]))
		resHets = res@gatk[training.site == T,]
		resHets$Sample = samps[x]
		resHets = resHets[sample(.N, rand_n)]
		gc()
		resHets = resHets[,c('af', 'Sample')]
		print(x)
		resHets
	})

	print('Combine hetSNPs across samples')
	cmnCols = Reduce(intersect, lapply(tmpL, function(x){colnames(x)}))
	alldf = lapply(tmpL, function(x){x[, ..cmnCols]}) %>% do.call(rbind, .)

	plt = gghistogram(alldf, x = "af",
	   add = "mean",
	   fill = "Sample", rug = F,
	   add_density = TRUE) +
	xlab('hetSNP_AF') + ylab('Count') +
	theme(axis.text.x = element_text(size=20),
		axis.text.y = element_text(size=20), legend.position = 'right',
		axis.title = element_text(size=20)) +
		rotate_x_text(90)

	print('Plotting allelic balance at hetSNPs')
	pdf(paste0(prefix, '_hetSNPs_AF_SCAN2.pdf'))
	print( plt )
	dev.off()

	return(plt)
}

plot_indel_mutprof = function(mutdf, genome_build, ref_genome,
				prefix = '', hg = 7, wd = 20, add_chr = F,
				ylims = c(0,200), oldcols = c('chr', 'pos', 'refnt', 'altnt')){

	setnames(mutdf, old = oldcols, new = c("chr","start","REF","ALT"), skip_absent=TRUE)

	if(add_chr == T){mutdf$chr = paste0('chr', mutdf$chr)}

	mutdf$end = mutdf$start
	mutgr = makeGRangesFromDataFrame(mutdf, keep.extra.columns=T)
	GenomeInfoDb::genome(mutgr) = genome_build

	indel_grl = get_indel_context(mutgr, ref_genome)

	indel_counts = count_indel_contexts(indel_grl)

	plt = plot_indel_contexts(indel_counts, condensed = TRUE)

	pdf(paste(prefix, 'MutationSignatures_Indel.pdf', sep = '_'), width = wd, height = hg)
	print(
	plt + theme(axis.text.x = element_text(size=15,face = 'bold'),
			axis.text.y = element_text(size=15,face = 'bold'),
			axis.title = element_text(size=15,face = 'bold')) +
	theme(text=element_text(size=15), axis.text.y = element_text(face = 'bold')) +
	rotate_x_text(90) + ylim(ylims)
	)
	dev.off()
}

scan2_to_pass = function(scan2_files, samps, ncores = 5, type = 'rda',
			passColumn = 'pass', subsample = F, subsampleN){

	suppressPackageStartupMessages(require(dplyr))
	suppressPackageStartupMessages(require(parallel))

	passL = mclapply(1:length(scan2_files), mc.cores = ncores, function(x){
		
		if(type == 'rda'){
			res = get(load(scan2_files[x]))
		} else if(type == 'rds'){
			res = readRDS(scan2_files[x])
		}
		resMuts = res@gatk[get(passColumn) == TRUE, ]
		if(subsample == T){
			print(paste0('Subsampling... ', samps[x]))
			resMuts <- resMuts[sample(.N, min(subsampleN, .N))]
		}
		resMuts$Sample = samps[x]
		colnames(resMuts)[colnames(resMuts) == samps[x]] = 'sc.gt'
		print(paste0('Loaded cell ', x))
		resMuts})

	cmnCols = Reduce(intersect, lapply(passL, function(x){colnames(x)}))
	pass_df = lapply(passL, function(x){x[, ..cmnCols]}) %>% do.call(rbind, .)
	return(pass_df)
}


scan2_to_burden = function(scan2_files, samps, ncores = 5, type = 'rda'){

	suppressPackageStartupMessages(require(tidyverse))
	suppressPackageStartupMessages(require(dplyr))
	suppressPackageStartupMessages(require(parallel))

	newbL = mclapply(1:length(scan2_files), mc.cores = ncores, function(x){
		if(type == 'rda'){
			res = get(load(scan2_files[x]))
			gburden = rbind(res@mutburden$snv[2,] %>% add_column(muttype = 'snv'),
			res@mutburden$indel[2,] %>% add_column(muttype = 'indel'))
			gburden$Sample = samps[x]
			gc()
			print(paste0('Loaded cell ', x))
			gburden
		} else if(type == 'rds'){
			res = readRDS(scan2_files[x])
			gburden = rbind(res@mutburden$snv[2,] %>% add_column(muttype = 'snv'),
			res@mutburden$indel[2,] %>% add_column(muttype = 'indel'))
			gburden$Sample = samps[x]
			gc()
			print(paste0('Loaded cell ', x))
			gburden
		}
		
	})

	cmnCols = Reduce(intersect, lapply(newbL, function(x){colnames(x)}))
	gburden_df = lapply(newbL, function(x){x[, cmnCols]}) %>% do.call(rbind, .)
	return(gburden_df)
}

plot_mutprof_from_vcf = function(mutdf, genome_build, ref_genome,
				prefix = '', hg = 7, wd1 = 14, wd2 = 10,
				add_chr = F, ymax1 = 0.2, ymax2 = 1, plot = T,
				oldcols = c('chr', 'pos', 'refnt', 'altnt')){

	setnames(mutdf, old = oldcols, new = c("chr","start","REF","ALT"), skip_absent=TRUE)

	if(add_chr == T){mutdf$chr = paste0('chr', mutdf$chr)}

	mutdf$end = mutdf$start
	mutgr = makeGRangesFromDataFrame(mutdf, keep.extra.columns=T)
	GenomeInfoDb::genome(mutgr) = genome_build
	
	mut_mat <- mut_matrix(vcf_list = mutgr, ref_genome = ref_genome)

	if(plot == T){

		pdf(paste(prefix, 'MutationSignatures_TripleNucleotide.pdf', sep = '_'), width = wd1, height = hg)
		print( plot_96_profile2(mut_mat, plot_size_scale = 2, condensed = T, ymax = ymax1) )
		dev.off()
	}

	type_occurrences <- mut_type_occurrences(mutgr, ref_genome)
	
	# Combine two specific groups (e.g., "G1" and "G2")
	df = as.data.frame(t(type_occurrences))
	df$ctgroup = rownames(df)
	df <- df %>%
	  mutate(group_combined = ifelse(ctgroup %in% c("C>T at CpG", "C>T other"), "C>T", ctgroup))
	df = df[df$ctgroup != 'C>T',]
	df$relcount = df$My_sample / sum(df$My_sample)

	# Plot
	cols = c('skyblue', 'black', 'pink', 'red', 'gray', 'darkolivegreen3', 'mistyrose')

	if(plot == T){

		pdf(paste(prefix, 'MutationSignatures_SingleNucleotide.pdf', sep = '_'), width = wd2, height = hg)
		print(	

		ggplot(df, aes(x = group_combined, y = relcount, fill = ctgroup)) +
		geom_bar(stat = "identity", position = "stack") +
		theme_classic() +
		xlab('') + ylab('Relative Contribution') +
		scale_fill_manual(values = cols) +
		theme(axis.text.x = element_text(size=20),
			axis.text.y = element_text(size=20),
			axis.title = element_text(size=20)) +
		ylim(0, ymax2) +
		rotate_x_text(90)

		)
		dev.off()
	}

	return(mut_mat)
}


extract_basecall_and_read = function(bam, outfile, chr, pos){

	# Define your parameters
	bam_file <- bam
	chromosome <- chr
	position <- pos

	python_script <- "/home/ch253794/SCRIPTS/Extract_readname_and_basecall2.py"  # Use the full path to your Python script

	# Construct the command to run the Python script
	command <- sprintf("python %s %s %s %s %s", python_script, bam_file, outfile, chromosome, position)

	# Execute the command
	system(command)
	
	# Read and return as data frame
	outfile_in = safe_read_table(outfile, colClasses = "character")
	return(outfile_in)
	
}


# Change behavior of plot_96_profile in mutationalpatterns
plot_96_profile2 = function(mut_matrix, colors = NA, ymax = 0.2, condensed = FALSE, plot_size_scale = 1) 
{
    freq <- full_context <- substitution <- context <- NULL
    if (is.na(colors)) {
        colors <- get("COLORS6", envir = asNamespace("MutationalPatterns"))
    }
    if (length(colors) != 6) {
        stop("Provide colors vector with length 6", call. = FALSE)
    }
    norm_mut_matrix <- apply(mut_matrix, 2, function(x) x/sum(x))
    tb <- norm_mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>% 
        dplyr::mutate(substitution = stringr::str_replace(full_context, 
            "\\w\\[(.*)\\]\\w", "\\1"), context = stringr::str_replace(full_context, 
            "\\[.*\\]", "\\.")) %>% dplyr::select(-full_context) %>% 
        tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", 
            values_to = "freq") %>% dplyr::mutate(sample = factor(sample, 
        levels = unique(sample)))
    if (condensed == TRUE) {
        width <- 1
        spacing <- 0
    }
    else {
        width <- 0.6
        spacing <- 0.5
    }
    plot <- ggplot(data = tb, aes(x = context, y = freq, fill = substitution, 
        width = width)) + geom_bar(stat = "identity", colour = "black", 
        size = 0.2) + scale_fill_manual(values = colors) + facet_grid(sample ~ 
        substitution) + ylab("Relative contribution") + coord_cartesian(ylim = c(0, 
        ymax)) + scale_y_continuous(breaks = seq(0, ymax, 0.1)) + 
        guides(fill = "none") + theme_bw() + theme(axis.title.y = element_text(size = 12*plot_size_scale, 
        vjust = 1), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12*plot_size_scale), 
        axis.text.x = element_text(size = 5*plot_size_scale, angle = 90, vjust = 0.5), 
        strip.text.x = element_text(size = 9*plot_size_scale), strip.text.y = element_text(size = 9*plot_size_scale), 
        panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing, 
            "lines"))
    return(plot)
}

