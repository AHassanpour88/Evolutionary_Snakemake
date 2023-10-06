#!/bin/bash

# set parameters
INSTALLCONDA="TRUE"
INSTALLMAMBA="TRUE"
ENVPATH="config/environment.yaml"
RANDOMFIELDS="https://github.com/tpook92/MoBPS/blob/master/RandomFieldsUtils_1.0.6.tar.gz"
MIRACULIX="https://github.com/tpook92/MoBPS/blob/master/miraculix_1.0.5.tar.gz"
MOBPSMAPS="https://github.com/tpook92/MoBPS/blob/master/MoBPSmaps_0.1.13.tar.gz"

# download and install latest miniconda if wanted
if [ ${INSTALLCONDA} == "TRUE" ]
then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc
fi

# install mamba
if [ ${INSTALLMAMBA} == "TRUE" ]
then
    conda install -y -c conda-forge mamba
fi

# create mobpsopti env and install MoBPS packages
mamba env create -f=${ENVPATH}
conda activate mobpsopti
R --vanilla -e "install.packages('${RANDOMFIELDS}',repos=NULL)"
R --vanilla -e "install.packages('${MIRACULIX}',repos=NULL)"
R --vanilla -e "remotes::install_github('https://github.com/tpook92/MoBPS', subdir = 'pkg')"
R --vanilla -e "install.packages('${MOPSMAPS}',repos=NULL)"

rm -f Miniconda3-latest-Linux-x86_64.sh
