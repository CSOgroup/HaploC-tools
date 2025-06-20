#!/bin/bash --login
eval "$(conda shell.bash hook)"

module load gcc miniconda3
conda activate nHapCUT2
export LC_ALL=C

##############################

while getopts ":d:k:" opt; do
  case $opt in
    d) wk_dir="$OPTARG"
    ;;
    k) k="$OPTARG"
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
        data_num=${SLURM_ARRAY_TASK_ID}
        Rscript $main_R $wk_dir $k $data_num
        ;;
    "bwa")
        n_cores=${SLURM_CPUS_PER_TASK}
        chunk_num=${SLURM_ARRAY_TASK_ID}
        Rscript $main_R $wk_dir $k $chunk_num
        ;;
    "hic")
        chunk_num=${SLURM_ARRAY_TASK_ID}
        Rscript $main_R $wk_dir $k $chunk_num
        ;;
    "rep")
        chunk_num=${SLURM_ARRAY_TASK_ID}
        $repair_sh -d $wk_dir -c $chunk_num
        ;;
    "sec")
        mkdir -p $wk_dir/mega/aligned/
        chunk_num=${SLURM_ARRAY_TASK_ID}
        Rscript $main_R $wk_dir $k $chunk_num
        ;;
    "intg")
        conda activate HapCUT2
        chr_num=${SLURM_ARRAY_TASK_ID}
        Rscript $shapeit_R $wk_dir $chr_num
        for maxIS_Mb in ${maxIS[@]}; do Rscript $hapcut_R $wk_dir $chr_num $maxIS_Mb & done; wait
        ;;
    "distr" | "opt" | "mates" | "snp" | "mega")
        Rscript $main_R $wk_dir $k
        ;;
    "htrans")
        Rscript $htrans_R $wk_dir
        ;;
    "htransplot")
        Rscript $htrans_plot_R $wk_dir
        ;;
    *)
        echo Invalid value of k: $k
        # Add any additional handling for an invalid value of k, if needed.
        ;;
esac
