#!/bin/bash --login
export LC_ALL=C

##############################

usage() {
    echo "Usage: haploc_cpu.sh [-d wk_dir] [-h]"
    echo "  -d    Set the working directory (wk_dir)."
    echo "  -x    Also include downstreams analysis (diffIns, diffComp and HaploCNV)."    
    echo "  -h    Display this help and exit."
}

for arg in "$@"; do [[ $arg == "--help" ]] && usage && exit 0; done
downstream=false

##############################

while getopts ":d:h:x:" opt; do
  case $opt in
    d) wk_dir="$OPTARG"
    ;;
    x) downstream="$OPTARG"
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

##############################

config_f=${wk_dir}/config/config.txt
ID=$(basename ${wk_dir}); echo $ID
sizeGb=$(grep '^sizeGb=' $config_f | cut -d '=' -f 2)
thread4bwa=$(grep '^thread4bwa=' $config_f | cut -d '=' -f 2)
repo_dir=$(grep '^repo_dir=' $config_f | cut -d '=' -f 2)
HaploC_sh=$repo_dir/HaploC/main/HaploC_helper_cpu.sh
downstream_sh=$repo_dir/bin/downstreams.sh

##############################

mkdir -p ${wk_dir}/log_file
echo ${wk_dir}

###############################

run_HaploC()
{

	############################## STEP0: split fastq

	k=split
	log_f=${wk_dir}/log_file/log.$k
	data_source_f=${wk_dir}/log_file/data_source.txt
	n_data=$(find ${wk_dir} -type f -name '*1.fastq.gz' | wc -l)
	find ${wk_dir} -type f -name '*_1.fastq.gz' -exec basename {} \; | sed 's/_1.fastq.gz//' > $data_source_f
	for ((id=1; id<=n_data; id++)); do $HaploC_sh -d ${wk_dir} -k $k -m $id; done

	############################## STEP1: distribute data to chunks

	k=distr
	$HaploC_sh -d ${wk_dir} -k $k

	##############################

	n_chunks_f=$wk_dir/log_file/n_chunks.txt
	while [ ! -f $n_chunks_f ]; do sleep 10; done
	n_chunks=$(cat "$n_chunks_f")

	############################## STEP2: align and so on	

	k=bwa
	for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk & done; wait
		
	############################## STEP3: proceed to generate hic

	k=hic
	for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk & done; wait

	############################## STEP3x: generate mega / do not do mega

	k=mega
	$HaploC_sh -d ${wk_dir} -k $k

	############################## STEP4: call snps, do not depend on previous step

	k=snp
	$HaploC_sh -d ${wk_dir} -k $k

	############################## STEP5: intersect SNP and gzip

	k=sec
	for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk & done; wait

	############################### STEP6: repair

	k=rep
	for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk & done; wait

	############################### STEP7: intg

	k=intg
	for chr in $(seq 1 23); do $HaploC_sh -d ${wk_dir} -k $k -c $chr & done; wait

	############################### STEP8: opt

	k=opt
	$HaploC_sh -d ${wk_dir} -k $k

	############################### STEP9: mates

	k=mates
	$HaploC_sh -d ${wk_dir} -k $k

	###############################

	[ "$downstream" = "false" ] && exit

	k=diffIns
	bin_size=25000	
	$downstream_sh -d $wk_dir -k $k -s $bin_size

	k=diffComp
	$downstream_sh -d $wk_dir -k $k

	k=HaploCNV
	bin_size=100000
	$downstream_sh -d $wk_dir -k $k -s $bin_size
}

##############################

run_HaploC $ID &

##############################
