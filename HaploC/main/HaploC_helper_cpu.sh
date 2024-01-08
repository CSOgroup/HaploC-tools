#!/bin/bash --login
# set -euo pipefail

conda activate nHapCUT2
export LC_ALL=C

##############################

while getopts ":d:k:m:n:c:" opt; do
  case $opt in
    d) wk_dir="$OPTARG"
    ;;
    k) k="$OPTARG"
    ;;
    m) data_num="$OPTARG"
    ;;
    n) chunk_num="$OPTARG"
    ;; 
    c) chr_num="$OPTARG"
    ;;       
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

##############################

config_f=$wk_dir/config/config.txt
repo_dir=$(grep '^repo_dir=' $config_f | cut -d '=' -f 2)
HaploC_dir=${repo_dir}/HaploC/

main_R=${HaploC_dir}/main/HaploC.R
htrans_R=${HaploC_dir}/H-trans-evaluation/H_trans_nOpt.R
htrans_plot_R=${HaploC_dir}/H-trans-evaluation/H_trans_plot.R
shapeit_R=${HaploC_dir}/HapCUT2/intg_shapeit.R
hapcut_R=${HaploC_dir}/HapCUT2/intg_hapcut.R
repair_sh=${HaploC_dir}/HapCUT2/repair.sh
maxIS_tmp=$(grep '^maxIS=' $config_f | cut -d '=' -f 2)
IFS=',' read -r -a maxIS <<< "$maxIS_tmp"

##############################

case $k in
    "split")
        ID="$(basename -- $wk_dir)"
        Rscript $main_R $wk_dir $k $data_num
        ;;
    "bwa")
        Rscript $main_R $wk_dir $k $chunk_num
        ;;
    "hic" | "rep")
        Rscript $main_R $wk_dir $k $chunk_num
        ;;
    "snp")
        Rscript $main_R $wk_dir $k
        ;;
    "sec")
        mkdir -p $wk_dir/mega/aligned/
        Rscript $main_R $wk_dir $k $chunk_num
        ;;  
    "intg")
        conda activate HapCUT2
        Rscript $shapeit_R $wk_dir $chr_num
        for maxIS_Mb in ${maxIS[@]}; do Rscript $hapcut_R $wk_dir $chr_num $maxIS_Mb & done; wait 
        ;;
    "distr" | "opt" | "mates" | "mega")
        Rscript $main_R $wk_dir $k
        ;;
    *)
        echo Invalid value of k: $k
        # Add any additional handling for an invalid value of k, if needed.
        ;;
esac
