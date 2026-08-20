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
conda create -n scan2 -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2

# Once the environment is created, it must be activated every time
# you wish to use SCAN2-Zoo.
conda activate scan2-zoo
```



### Internal R reference genome
```
conda install -c conda-forge -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg38
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
SCAN2-Zoo relies on high converage (>15X) and a panel of cells from
different subjects for accurate somatic SNV and indel calling. We
provide a minimal example here for users to test whether the tool
runs with the expected output files. Somatic mutation calls of this
tutorial are not accurate. For a more detailed tutorial with
explanations, please follow our vignette: XYZ

## Download the BAM files
This is a glutamatergic neuron from a marmoset sample and its matched
bulk WGS (aligned to calJac4). Both are subset to ~2X genomic coverage
and limited to chromosome 20 for tutorial purposes.

wget XYZ

## Download the panel files
Since building a panel requires multiple samples, here we provide the
pre-built panel files to enable full run of the third step for tutorial
purposes.

wget XYZ

## Download the homopolymer regions
Refinement step removes somatic mutation calls that map to homopolymers
longer than 5 nucleotides. Here we provide this file for calJac4, please
see the vignette for generating homopolymer regions from any genome.

wget XYZ

## Generate a regions file
To parallelize genotyping with gatk, genome is split into a user defined
number of regions. Here we will use calJac4 chromosome 20 split into 3
intervals. This can be increased for larger datasets (typically 100-300
regions for the entire genome).

wget XYZ

# Run SCAN2-Zoo

## Genotyping

scan2-zoo \
  -d PTA_230 \
  --snakefile ${SCAN2_dir}/snakemake/Snakefile \
  config \
  --gatk gatk4_gvcf \
  --scripts ${SCAN2_dir}/scripts \
  --bulk-bam ${BULK} \
  --sc-bam ${SC1} \
  --genome 'calJac4' \
  --ref /home/ch253794/workdir/REFERENCES/BWA_REFERENCES/calJac4/calJac4.fa \
  --regions-file /home/ch253794/workdir/REFERENCES/GenomicIntervals/calJac4_intervals.txt



## Run SCAN2
Replace instances of /path/to/... with the paths
downloaded above. This demo runs in ~15 minutes on a 20 core
machine by restricting analysis to a 2 MB segment of chr22.
```
# Creates a directory for the analysis
scan2 -d demo init
cd demo

# Configure analysis parameters. Can be run multiple times to change
# parameters.
scan2 config \
    --verbose \
    --ref /path/to/human_g1k_v37_decoy.fasta \
    --dbsnp /path/to/dbsnp_138_b37.vcf \
    --shapeit-refpanel  /path/to/1000GP_Phase3 \
    --regions 22:30000000-30999999,22:31000000-31999999 \
    --bulk-bam /path/to/hunamp.chr22.bam \
    --sc-bam /path/to/il-11.chr22.bam 
    --sc-bam /path/to/il-12.chr22.bam 

# Attempts to validate the chosen parameters and input files. Not
# exhaustive, but helpful to prevent errors.
scan2 validate

# Run the analysis. Set the number of cores you wish you use; for
# more advanced parallelization on a cluster or the cloud, see the
# sections below.
scan2 run --joblimit <N cores>
```
See `scan2 -h` for more details on configuration options and runtime arguments.

After SCAN2 completes, single sample results are stored in 
Rdata format `demo/call_mutations/[single_cell_sample_name]/scan2_object.rda`.

Run R and load the RData file. IMPORTANT! Use the R installation in the SCAN2 conda
environment for the SCAN2 R library.
```
$ R   
R> library(scan2)
# The two single cell samples are 'h24' (il-11) and 'h25' (il-12)
R> load('demo/call_mutations/h25/scan2_object.rda')
```

We now provide an S4 class called 'SCAN2' that handles both the
genotyping logic and allows convenient user interaction. Each
single cell in a SCAN2 run will have a 'results' object saved in the
corresponding scan2_object.rda file.


The S4 class provides a summary of several features of the data.
```
R> results
# SCAN2 
#   Region:
GRanges object with 22 ranges and 0 metadata columns:
       seqnames      ranges strand
          <Rle>   <IRanges>  <Rle>
   [1]        1 1-249250621      *
   [2]        2 1-243199373      *
   [3]        3 1-198022430      *
   [4]        4 1-191154276      *
   [5]        5 1-180915260      *
   ...      ...         ...    ...
  [18]       18  1-78077248      *
  [19]       19  1-59128983      *
  [20]       20  1-63025520      *
  [21]       21  1-48129895      *
  [22]       22  1-51304566      *
  -------
  seqinfo: 22 sequences from GRCh37.p13 genome
#   Genome: hs37d5 
#   Single cell ID: h25 
#   Bulk ID: hunamp 
#   GATK: 2517 raw sites
#   AB model training hSNPs: 753 phasing: 0|1=411, 1|0=342
#       VAF correlation between neighboring hSNPs:
#           <10 bp 1 <100 bp 0.983 <1000 bp 0.897 <10 kbp 0.806 <100 kbp 0.203 
#         36 resampled hSNPs
#         7 resampled hIndels
#   AB model parameters:
#       average (over chromosomes): a=-2.374, b=669.314, c=1.125, d=10193.908
#   Allele balance:
#       mean (0 is neutral): -0.073 
#       uncertainty (Q25, median, Q75): 0.372 0.558 1.162 
#       mean at training hSNPs: 0.058 
#   Mutation models: computed
#   CIGAR data: single cell: 2517 sites, bulk: 2517 sites
#   Static filters: 
#          snv:       12 retained     1873 removed        0 NA
#        indel:        7 retained      625 removed        0 NA
#   FDR prior data: 
#          snv:       25 candidates      657 germline hets        4 max burden
#        indel:       13 candidates      108 germline hets        0 max burden
#   Depth profile: 
#          snv: genome <:   min. bulk DP 1.3%,   min. sc DP 8.6%,   either 9.4%
#        indel: genome <:   min. bulk DP 1.3%,   min. sc DP 15.4%,   either 15.9%
#   Somatic mutation calls: mode=new, target.fdr=0.010
#          snv:        0 called       11 resampled training calls
#        indel:        0 called        1 resampled training calls
#       ALL indel calls have been suppressed (insufficient number of single cells in cross-sample filter
#   Somatic mutation burden: 
#          snv:     NA somatic,   NA% sens,   NA callable Gbp,   NA muts/haploid Gbp,   NA muts per genome
#        indel:     NA somatic,   NA% sens,   NA callable Gbp,   NA muts/haploid Gbp,   NA muts per genome
#       indel mutation burden estimates ARE NOT VALID (cross-sample panel insufficiency)!
#   Mutation rescue by signature: (not computed)
```

This SCAN2 object contains only VAF-based somatic mutation calls. For
signature-based rescue, see `scan2 rescue -h`. Further, this SCAN2 object
will not contain indel calls because no cross sample filter was provided
(see `scan2 config --help`, `--cross-sample-panel`).

There should be no called somatic SNVs in this demo dataset.
```
R> snvs <- results@gatk[pass == TRUE & muttype == 'snv']
R> snvs
Empty data.table (0 rows and 73 cols): chr,pos,dbsnp,refnt,altnt,mq...
```


# Parallelization
In a practical setting, parallelization will be required. SCAN2 leverages
Snakemake to offer parallelization via any of the following:
* A single machine with multiple cores. To do this, increase the
  `--joblimit` parameter but do not invoke cluster or cloud
  arguments.  A total memory limit can also be set via
  `--memlimit`, which should be supplied in megabytes.
* Clusters with distributed resource management software (e.g., SLURM,
  LSF, GridEngine). If your resource management software supports
  DRMAA (https://www.drmaa.org), the authors recommend using the
  `--drmaa`. **However, additional libraries may be necessary to
  interface with the DRMAA wrapper.**
    * E.g., in the SCAN2 publication, a SLURM cluster was accessed via
      the slurm-drmaa 1.1.1 package (https://github.com/natefoo/slurm-drmaa). This
      package was maintained by local cluster system admins.
      An example command line for running SCAN2 with DRMAA follows. Note the `--mem`
      argument, which uses Snakemake's `{resources.mem}` placeholder to allow each
      individual job in the SCAN2 pipeline to specify the amount of memory needed and
      the `{threads}` placeholder.
        ```
        scan2 run --joblimit 100 --drmaa ' -p YOUR_QUEUE_NAME -A YOUR_ACCOUNT -t 10:00 -c {threads} --mem={resources.mem}'
        ```
    * If your scheduler is not DRMAA-compatible (or if the appropriate
      DRMAA interface is unavailable), Snakemake's `--cluster` option
      offers similar functionality to `--drmaa`, but with fewer features.
* Cloud environments.

See Snakemake's documentation for more details on cluster and cloud execution: https://snakemake.readthedocs.io/en/stable/.

