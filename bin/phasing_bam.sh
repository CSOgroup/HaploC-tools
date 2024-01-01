#!/bin/bash --login
export LC_ALL=C
conda activate nHapCUT2

##############################

usage() {
    echo "Usage: haploc_cpu.sh [-d wk_dir] [-h]"
    echo "  -d    Set the working directory (wk_dir)."
    echo "  -h    Display this help and exit."
}

for arg in "$@"; do [[ $arg == "--help" ]] && usage && exit 0; done

##############################

while getopts ":d:b:h:" opt; do
  case $opt in
    d) wk_dir="$OPTARG"
    ;;
    b) bam2phase="$OPTARG"
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

###################################

config_f=$wk_dir/config/config.txt
repo_dir=$(grep '^repo_dir=' $config_f | cut -d '=' -f 2)
juicerDir=$repo_dir/HaploC/dependancies/juicer

###################################

bam_base=$(basename "$bam2phase" .bam)
bam_dir=$(dirname "$bam2phase")   

echo $bam_base

###################################

# bin_size=10000
# sbin_size=$((bin_size * 100))
# unphased_bw=$bam_dir/${bam_base}.unphased.bw

# echo $unphased_bw
# if [ ! -f "$unphased_bw" ]; then
	# samtools index -@ 20 $bam2phase
	# bamCoverage --bam $bam2phase -o $unphased_bw --binSize 1000 --smoothLength 10000 -p 50 --outFileFormat bigwig --normalizeUsing CPM ## no need for small binsize for unphased
# fi

###################################

phasing_bam_fun()
{
  	chr_num=$1

	chr_region=chr$chr_num
	echo chr$chr_num

	sam_region=$bam_dir/${bam_base}_${chr_region}.delete
	output_dip=${sam_region}.dip

	samtools view $bam2phase -@ 2 $chr_region > $sam_region
	samtools view -H $bam2phase > ${sam_region}.mat.sam ## header
	samtools view -H $bam2phase > ${sam_region}.pat.sam ## header

	phased_bedgraph_mat=$bam_dir/$bam_base.chr${chr_num}.mat.bedgraph
	phased_bedgraph_pat=$bam_dir/$bam_base.chr${chr_num}.pat.bedgraph
	phased_bedgraph_mat_minus_pat=$bam_dir/$bam_base.chr${chr_num}.diff.bedgraph

	chr_pos_f=$wk_dir/mega/HapCut/phased_hap_intg/hic_chr${chr_num}_opt.chr_pos
	snp_mates_f=$wk_dir/mega/HapCut/phased_hap_intg/hic_chr${chr_num}_opt.pat_mat

    if [ ! -e "$chr_pos_f" ]; then
        exit 0
    fi

	cat $sam_region | $juicerDir/scripts/common/diploid_on_sam_YL.pl -s $chr_pos_f -o $snp_mates_f > ${output_dip}
	grep 'maternal' ${output_dip} | cut -d" " -f2- >> ${sam_region}.mat.sam
	grep 'paternal' ${output_dip} | cut -d" " -f2- >> ${sam_region}.pat.sam

	samtools view ${sam_region}.pat.sam -b > ${sam_region}.pat.bam
	samtools view ${sam_region}.mat.sam -b > ${sam_region}.mat.bam
}

###################################

chr_nums=(X $(seq 1 1 22))

# if [ ! -f "$bam_dir/${bam_base}.mat.bam" ]; then
	for chr_num in ${chr_nums[@]}; ## this can be parallelized
	do
		phasing_bam_fun $chr_num &
	done
	wait

	###################################

	samtools merge -f $bam_dir/$bam_base.pat.bam $bam_dir/$bam_base*delete.pat.bam
	samtools index $bam_dir/${bam_base}.pat.bam

	samtools merge -f $bam_dir/$bam_base.mat.bam $bam_dir/$bam_base*delete.mat.bam
	samtools index $bam_dir/${bam_base}.mat.bam

	rm $bam_dir/${bam_base}*delete*
# fi

###################################

# bin_size=1000
# bamCoverage --bam $bam_dir/$bam_base.mat.bam -o $bam_dir/$bam_base.mat.bedgraph --binSize $bin_size  -p 30 --outFileFormat bedgraph
# bamCoverage --bam $bam_dir/$bam_base.pat.bam -o $bam_dir/$bam_base.pat.bedgraph --binSize $bin_size  -p 30 --outFileFormat bedgraph

# exit

# bamCoverage --bam $bam_dir/$bam_base.mat.bam -o $bam_dir/$bam_base.mat.bw --binSize $bin_size --smoothLength ${sbin_size} -p 30 --outFileFormat bigwig
# bamCoverage --bam $bam_dir/$bam_base.pat.bam -o $bam_dir/$bam_base.pat.bw --binSize $bin_size --smoothLength ${sbin_size} -p 30 --outFileFormat bigwig

# bigwigCompare -b2 $bam_dir/$bam_base.mat.bw -b1 $bam_dir/$bam_base.pat.bw -bs $bin_size -p 30 --outFileFormat bigwig -o $bam_dir/$bam_base.pat_minus_mat.bw --operation subtract
# bigwigCompare -b2 $bam_dir/$bam_base.mat.bw -b1 $bam_dir/$bam_base.pat.bw -bs $bin_size -p 30 --outFileFormat bigwig -o $bam_dir/$bam_base.pat_minus_mat.log2.bw

# bigwigCompare -b2 $bam_dir/$bam_base.mat.bw -b1 $bam_dir/$bam_base.pat.bw -bs $bin_size -p 30 --outFileFormat bedgraph -o $bam_dir/$bam_base.pat_minus_mat.bedgraph --operation subtract
# bigwigCompare -b2 $bam_dir/$bam_base.mat.bw -b1 $bam_dir/$bam_base.pat.bw -bs $bin_size -p 30 --outFileFormat bedgraph -o $bam_dir/$bam_base.pat_minus_mat.log2.bedgraph
