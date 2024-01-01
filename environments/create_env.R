########################### Enviroment for HapCUT2

conda env remove -n HapCUT2
conda create -y -n HapCUT2 python=2.7 ## need to install python 2.7
conda activate HapCUT2
conda install -y -c bioconda hapcut2
conda install -y -c r r-base

# packages = c('doParallel')
# for (package in packages) {
#   if (!requireNamespace(package, quietly = TRUE)) {
#     print(package)
#     try(install.packages(package, dependencies = TRUE))
#   }
# }

########################### Enviroment for n_phasing

## python and R

conda env remove -n nHapCUT2
## These tools need be already installed (version not required)
conda create -y -n nHapCUT2
conda activate nHapCUT2
conda install -y -c bioconda samtools=1.19
conda install -y -c bioconda bwa
conda install -y -c bioconda bamtools
conda install -y -c bioconda vcftools
conda install -y -c bioconda vcflib=1.0.0_rc2
conda install -y -c bioconda bcftools
conda install -y -c bioconda bbmap
conda install -y -c r r-base ## no need to install the latest version (4.3.xx), from 3.6.xx will be fine (not fully tested)

## R packages
# data.table, doParallel, GenomicRanges, ggplot2, gplots, grid, Matrix, tools

########################### Enviroment for downstream analsis

## seg fault of R: https://stackoverflow.com/questions/76615222/r-crashes-immediately-on-startup-caught-segfault-address-nil-cause-memory
su -l ${USER}
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GenomicRanges")

packages = c('mgcv', 'rpart', 'R.utils', 'doParallel', 'ape', 'dendextend', 'fitdistrplus', 'igraph', 'Matrix', 'rARPACK', 'data.table', 'fields', 'GenomicRanges', 'strawr', 'ggplot2', 'Rcpp', 'RcppArmadillo')
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    print(package)
    try(install.packages(package, dependencies = TRUE))
  }
}

# install.packages(dependencies = TRUE, "/mnt/ptemp/Yuanlong/1.Projects/15.Phasing/0.Scripts/CALDER2/", repos = NULL, type = "source")
install.packages(dependencies = TRUE, "/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/HaploC-tools/CALDER2/", repos = NULL, type = "source")

# find ./ -type f -name "*.R" -exec grep -Ril 'chrom_sizes = ' {} \;
# find ./ -type f -exec sed -i 's|HaploC/header_v3.R|header/HaploC_header.R|g' {} +
# find ./ -type f -exec sed -i 's|n_chunks.txt|log_file/n_chunks.txt|g' {} +
# find ./ -type f -name "*.R" -exec grep -Ril 'chrom_sizes = ' {} \;
# find ./ -type f -exec sed -i 's|forceSymmetric|my_forceSymmetric|g' {} +
# find ./ -type f -name "*.R" -exec grep -Ril 'job_ids.txt' {} \;

# scp -r yliu5@curnagl.dcsr.unil.ch://work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/ /mnt/ptemp/Yuanlong/1.Projects/15.Phasing
# scp -r yliu5@curnagl.dcsr.unil.ch://scratch/yliu5/1.Projects/15.Phasing/1.Data/small_hic /mnt/ptemp/Yuanlong/1.Projects/15.Phasing

find ./ -type f -name "*.sh" -exec sed -i 's|HapCUT2|HapCUT2|g' {} +
find ./ -type f -exec sed -i 's|HaploC|HaploC|g' {} +

find ./ -type f -name "*.sh" -exec sed -i 's|n_bwa_thread|thread4bwa|g' {} +



########################### HapCUT2

conda activate HapCUT2
conda env export > env_HapCUT2.yml

conda activate nHapCUT2
conda env export > env_nHapCUT2.yml


export ZENODO_TOKEN=aGClYdvJn6ihjQ5FjLG8gN4lI5C8upNec4UPOhPgsTAHo00YneYWpTZ1gqgG
conda activate intg
conda install -c conda-forge jq
zenodo-upload/zenodo_upload.sh 10446020 genomicData.tar.gz 
zenodo-upload/zenodo_upload.sh 10446020 /scratch/yliu5/1.Projects/15.Phasing/1.Data/demo_data.tar.gz

tar -czvf demo_data.tar.gz demo_data
