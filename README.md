# SCAN2-Zoo
Genotyper for **somatic** SNV and indel discovery in PTA-amplified single cells
from any diploid organism.

## Important Notes

SCAN2-Zoo should only be applied to diploid autosomes (i.e., human autosomes 1-22).
Sex chromosomes are not supported. SCAN2-Zoo has so far been tested on the following
organisms: human, chimpanzee, rhesus macaque, marmoset, ferret and hybrid mouse
(BALB/c X C57BL/6).

SCAN2-Zoo aggressively removes mutations with any supporting read in the matched
bulk WGS. The idea behind this is to remove any germline variants but this will
also remove early developmental somatic mutations, especially if the matched bulk
WGS has high depth (>60X) and is from the same tissue as the single-cell. Therefore
we do not recommend SCAN2-Zoo for early developmental somatic mutation calling.

SCAN2-Zoo is essentially a modification of SCAN2 that was developed to call
mutations from human cells. We refer the readers to the original tool
if they are only interested in calling mutations from human cells:
<https://github.com/parklab/SCAN2>


# License
SCAN2-Zoo is freely available for non-commercial use.

# Installation
SCAN2-Zoo is distributed as a conda package. Installation requires the conda
package management tool and a Linux-flavored OS.

## Install SCAN2-Zoo
Create a conda environment for SCAN2-Zoo and install necessary packages
```
conda build -c conda-forge -c bioconda -c ecaglayan recipe/scan2-zoo

# Once the environment is created, it must be activated every time
# you wish to use SCAN2-Zoo.
conda activate scan2-zoo
```

# Overview of SCAN2-Zoo
SCAN2-Zoo runs in three steps; genotyping, panel building and refinement.

Genotyping: Genotypes all single-cell and bulk WGS from the same subject
and gives a scan2 object for each cell that has the genotyping statistics
for both the single-cell and its matched bulk (object@gatk). This step
also builds an allelic balance curve based on the maximum VAF of heterozygous
SNPs and creates statistics to be used in somatic mutation calling. Although
this step will also call somatic SNVs, the results are not final and
should not be used.

Panel building: Combines genotypes of all single-cell and bulk WGS in the
entire dataset (i.e all subjects) and gives statistics on all variants to
be used in refinement.

Refinement: Uses additional filters and panel information for final SNV and
indel calls.

# Minimal SCAN2-Zoo run
SCAN2-Zoo relies on high coverage (>10X) and a panel of cells from
different subjects for accurate somatic SNV and indel calling. We
provide a minimal example here for users to test whether the tool
runs with the expected output files. Somatic mutation calls of this
tutorial are not final. All files are hosted in Zenodo:

ZENODO LINK

## The BAM files
This is a glutamatergic neuron from a marmoset sample and its matched
bulk WGS (aligned to calJac4). Both are either limited to chromosome 20 (for step 1) or
subset to ~10X whole genome coverage (for step 3) for tutorial purposes.

## The panel files
Since building a panel requires multiple samples, here we provide the
pre-built panel files of marmosets to enable full run for tutorial
purposes. These files are subset to chr20 for step 1.

## The homopolymer regions
Refinement step removes somatic mutation calls that map to homopolymers
longer than 5 nucleotides. Here we provide this file for calJac4 chr20.
Please see the vignette for generating homopolymer regions from any genome.

## The regions file
To parallelize genotyping with gatk, genome is split into a user defined
number of regions. Here we will use calJac4 chromosome 20 split into 5
intervals (for step 1). This can be increased for larger datasets (typically 100-300
regions for the entire genome).

## Download the marmoset genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/calJac4/bigZips/calJac4.fa.gz
gunzip calJac4.fa.gz

# Run SCAN2-Zoo

## Step 1: Genotype


```

# Change to your work directory that has sufficient storage

# Create a folder that for the outputs
mkdir -p scan2_zoo_mini

# Initialize SCAN2-Zoo for this run
scan2-zoo \
-d scan2_zoo_mini \
init

# Configure SCAN2-Zoo for this run
scan2-zoo \
-d scan2_zoo_mini \
config \
--gatk gatk3_gvcf \
--bulk-bam datasets/bulk_230_chr20_subset.bam \
--sc-bam datasets/sc_D10_230_chr20_subset.bam \
--genome 'calJac4' \
--ref /path/to/calJac4.fa \
--regions-file datasets/calJac4_intervals_chr20.txt

# Validate
cd scan2_zoo_mini
scan2-zoo validate

# Run
scan2-zoo run --joblimit 8


```


## Step 2: Build the cross-sample panel

We will skip panel building for this example and use the prebuilt panel.
Please refer to the detailed wiki for building a panel from multiple
samples.


## Step 3: Refine


```

scan2-zoo \
-d scan2_zoo_mini_refine \
init

# Add the scan2 object paths to a file. Normally this would be the files from step1 but
# in this case we will not use it since it is not genomewide. We have uploaded
# the genomewide calls under 'datasets' in Zenodo.

realpath datasets/scan2_object.rda > scan2_zoo_mini_refine/scan2_objects.txt

scan2-zoo \
	-d scan2_zoo_mini_refine \
	config \
	--verbose \
	--analysis refine \
	--scan2-object-file scan2_zoo_mini_refine/scan2_objects.txt \
	--indel-fraction-delta 0.5 \
	--genome 'calJac4' \
	--bsgenome 'BSgenome.Cjacchus.UCSC.calJac4' \
	--ref /path/to/calJac4.fa \
	--makepanel-metadata datasets/panel_marmoset.csv \
	--panel datasets/merge_mmq60.tab.gz \
	--homopolymers datasets/calJac4_homopolymers.bed \
	--filter-homopolymers \
	--target-fdr 0.01 \
	--panel-stats datasets/panel_230_panel.tab.gz \
	--bulk-bam datasets/bulk_230_subset.bam


cd scan2_zoo_mini_refine
scan2-zoo refine --snakemake-args "--configfile scan.yaml --config config_path=scan.yaml"


```


## Examine results

Final step of SCAN2-Zoo creates a file that has the mutation burden, mutations that pass
the filters and mutations that failed the final filters (in the 'refine' step) and which
filter they failed. SCAN2-Zoo also creates a SCAN2 object that contains statistics
for every variant in addition to many other details (scan2_object_refine.rda).
Here we will load the final mutation calls and mutation burden estimates as well as the
final SCAN2 object.

Run R and load the output files as indicated below. Make sure to use the R installation
in the SCAN2-Zoo conda environment for the SCAN2-Zoo R library.

```
$ R   
R> library(scan2.zoo)
R> res = readRDS('objects/230-D10/230-D10_final_burdens_and_calls.rds')
R> names(res)

R> # Get genomewide burden estimates
R> burden = res[['Burden']]
R> # The most accurate burden estimate is the second quantile (Q2) that
R> # corresponds to somatic calls sequenced at normal depth (between 25th and 75th percentile).
R> snv_burden = burden[burden$muttype == 'snv' & burden$quantile == 'Q2', 'burden']
R> print(snv_burden)
R> indel_burden = burden[burden$muttype == 'indel' & burden$quantile == 'Q2', 'burden']
R> print(indel_burden)

R> # Examine the final mutation calls
R> head(res[['Passes_snv']])
R> head(res[['Passes_indel']])

R> # The mutations filtered in the last 'refine' step can also be examined.
R> # column 'reason' indicates why it is filtered
R> head(res[['Filtered_snv']])
R> head(res[['Passes_indel']])

R> # Load the final SCAN2 object if you would like to check statistics on all
R> # variants called by GATK.
R> obj = get(load('objects/230-D10/230-D10_scan2_object_refine.rda'))
R> head(obj@gatk)


```





