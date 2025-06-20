#!/bin/bash --login
eval "$(conda shell.bash hook)"

export LC_ALL=C

##############################

while getopts ":d:c:" opt; do
  case $opt in
    d) wk_dir="$OPTARG"
    ;;
    c) chunk_num="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

#########################################################

config_f=$wk_dir/config/config.txt
genome_version=$(grep '^genome_version=' $config_f | cut -d '=' -f 2)
repo_dir=$(grep '^repo_dir=' $config_f | cut -d '=' -f 2)
HaploC_dir=$repo_dir/HaploC
HiC_repair_py=$HaploC_dir/HapCUT2/my_hic_repair.py

#########################################################

if [ $genome_version == mm10 ]; then
 chr_nums=({1..19} X)
fi

if [ $genome_version == hg19 ]; then
 chr_nums=({1..22} X)
fi

#########################################################

wk_dir_chunk=${wk_dir}/chunk${chunk_num}/
tmp_dir=${wk_dir_chunk}/hapcut_tmp/
hic_split_dir=${wk_dir_chunk}/hic_separated/

#########################################################

mkdir -p $hic_split_dir
mkdir -p $tmp_dir

#########################################################

input_mate1=$(find "${wk_dir_chunk}/splits/" -type f -name "*_R1.n_sorted.bam" | head -n 1)
input_mate2=$(find "${wk_dir_chunk}/splits/" -type f -name "*_R2.n_sorted.bam" | head -n 1)
hic_repaired=${wk_dir_chunk}/splits/hic_repaired.sam
hic_fixmate_sorted=${wk_dir_chunk}/splits/hic_fixmate_sorted.bam

######################################################### repair chimeric, use conda activate HapCUT2

# repair chimeric Hi-C read pairs (from off-center Hi-C linker)
conda run -n HapCUT2 python2.7 ${HiC_repair_py} ${input_mate1} ${input_mate2} ${hic_repaired} ${HaploC_dir}/dependancies/HapCUT2/

#########################################################

conda activate nHapCUT2
samtools fixmate -@ 10 ${hic_repaired} - | samtools sort -@ 10 -T ${tmp_dir} -o ${hic_fixmate_sorted} - && rm ${hic_repaired}
samtools index -@ 10 ${hic_fixmate_sorted}

######################################################### split into chrs

for chr_num in "${chr_nums[@]}";
do
    chr_char=chr$chr_num
    samtools view -b ${hic_fixmate_sorted} ${chr_char} -@ 10 -o ${hic_split_dir}/hic_${chr_char}.bam
    echo "chr$chr_num split done"
done

#########################################################

# fz = file.size( sprintf('%s/hic_%s.bam', hic_split_dir, chr_char) )/1024/1024 ## Mb
# # if(fz > 1)
# {
#     unlink(input_mate1)
#     unlink(input_mate2)
# }
