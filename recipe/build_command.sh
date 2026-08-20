#!/bin/bash

#conda-build -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2


#conda create -n build_env conda-build
#conda activate build_env

#conda install conda-build


conda build -c conda-forge -c bioconda -c ecaglayan recipe/scan2-zoo

