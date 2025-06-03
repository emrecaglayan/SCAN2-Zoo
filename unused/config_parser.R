#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)


if (!(length(args) %in% 2:3)){
    stop("usage: config_parser.R prefix")
}

prefix <- args[1]

# Step 1: List matching top-level PTA* directories
pref_dirs <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)
pref_dirs <- grep(paste0(prefix, ".*"), pref_dirs, value = TRUE)

# Step 2: Search for 'scan.yaml' directly inside each
yaml_files <- unlist(lapply(pref_dirs, function(d) {
  list.files(path = d, pattern = "^scan.yaml$", full.names = TRUE)
}))

# Step 3: Loop throught the config files and generate a data frame of bulk bam and scan2 objects
dfL = list()
for(i in 1:length(yaml_files)){

	config = yaml::read_yaml(yaml_files[i])
	bulkbam = config$bulk_bam
	scan2_objs = list.files(pattern = 'scan2_object_update.rda', path = basename(dirname(yaml_files[i])), recursive = T)
	scan2_objs = paste0(getwd(), '/', basename(dirname(yaml_files[i])), '/', scan2_objs)
	dfL[[i]] = data.frame(bulkbam = bulkbam, scan2_obj = scan2_objs)
}

df = do.call(rbind, dfL)


