#!/bin/bash
set -euo pipefail

# Executables
mkdir -p $PREFIX/bin
cp bin/scan2-zoo $PREFIX/bin/scan2-zoo
cp bin/joint_panel $PREFIX/bin/joint-panel

# Supporting R/Python/bash scripts
mkdir -p $PREFIX/lib/scan2-zoo
cp scripts/*.R scripts/*.py scripts/*.sh $PREFIX/lib/scan2-zoo/

# Snakemake pipeline files
cp snakemake/Snakefile $PREFIX/lib/scan2-zoo/
cp snakemake/snakefile.* $PREFIX/lib/scan2-zoo/

# Snakemake's own scripts subfolder
mkdir -p $PREFIX/lib/scan2-zoo/scripts
cp snakemake/scripts/*.R $PREFIX/lib/scan2-zoo/scripts/
