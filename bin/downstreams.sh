#!/bin/bash --login
export LC_ALL=C
conda activate nHapCUT2

##############################

usage() {
    echo "Usage: downstreams.sh [-d wk_dir] [-k module] [-s bin_size] [-h]"
    echo "  -d    Set the working directory (wk_dir)."
    echo "  -k    Set the module to run, should be one of: diffIns, diffComp, HaploCNV"
    echo "  -s    Set the bin size. Default value: 25000 for diffIns analysis, 100000 for HaploCNV analysis"
    echo "  -h    Display this help and exit."
}

for arg in "$@"; do [[ $arg == "--help" ]] && usage && exit 0; done

##############################

while getopts ":d:k:s:h" opt; do
  case $opt in
    d) wk_dir="$OPTARG"
    ;;
    k) k="$OPTARG"
    ;;
    s) bin_size="$OPTARG"
    ;;
    h) usage
       exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
       usage
       exit 1
    ;;
  esac
done

if [ -z $bin_size ]; then
    if [ "$k" = "diffIns" ]; then
        bin_size=25000
    elif [ "$k" = "HaploCNV" ]; then
        bin_size=100000
    fi
fi

##############################

config_f=$wk_dir/config/config.txt
ID=$(basename "$wk_dir"); echo $ID
sizeGb=$(grep '^sizeGb=' $config_f | cut -d '=' -f 2)
thread4bwa=$(grep '^thread4bwa=' $config_f | cut -d '=' -f 2)
repo_dir=$(grep '^repo_dir=' $config_f | cut -d '=' -f 2)

diffIns_R=$repo_dir/diffIns/insulation.R
HaploCNV_R=$repo_dir/HaploCNV/HaploCNV_main.R
diffComp_R=$repo_dir/diffComp/diffComp_main.R

##############################

case $k in
    "diffIns")
        Rscript $diffIns_R $wk_dir $bin_size
        ;; 
    "diffComp")
        Rscript $diffComp_R $wk_dir
        ;;       
     "HaploCNV")
        Rscript $HaploCNV_R $wk_dir $bin_size
        ;; 
    *)
        echo Invalid value of k: $k
        ;;
esac
