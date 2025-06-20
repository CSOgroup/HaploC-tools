#!/bin/bash --login
eval "$(conda shell.bash hook)"
export LC_ALL=C

##############################

usage() {
    echo "Usage: haploc_cpu.sh [-d wk_dir] [-x] [-h]"
    echo "  -d    Set the working directory (wk_dir)."
    echo "  -x    Include downstream analysis (diffIns, diffComp, and HaploCNV) when set as "true". Optional."
    echo "  -h    Display this help and exit."
}

for arg in "$@"; do [[ $arg == "--help" ]] && usage && exit 0; done
downstream=false

##############################

while getopts ":d:x:h" opt; do
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
    log_f=${wk_dir}/log_file/log.$k.txt
    data_source_f=${wk_dir}/log_file/data_source.txt
    n_data=$(find ${wk_dir} -type f -name '*1.fastq.gz' | wc -l)
    find ${wk_dir} -type f -name '*_1.fastq.gz' -exec basename {} \; | sed 's/_1.fastq.gz//' > $data_source_f
    for ((id=1; id<=n_data; id++)); do $HaploC_sh -d ${wk_dir} -k $k -m $id >> $log_f 2>&1; done

    ############################## STEP1: distribute data to chunks

    k=distr
    log_f=$wk_dir/log_file/log.$k.txt
    $HaploC_sh -d ${wk_dir} -k $k >> $log_f 2>&1

    ##############################

    n_chunks_f=$wk_dir/log_file/n_chunks.txt
    while [ ! -f $n_chunks_f ]; do sleep 10; done
    n_chunks=$(cat "$n_chunks_f")

    ############################## STEP2: align and so on

    k=bwa
    log_f=$wk_dir/log_file/log.$k.txt
    for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk >> $log_f 2>&1 & done; wait

    ############################## STEP3: proceed to generate hic

    k=hic
    log_f=$wk_dir/log_file/log.$k.txt
    for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk >> $log_f 2>&1 & done; wait

    ############################## STEP3x: generate mega / do not do mega

    k=mega
    log_f=$wk_dir/log_file/log.$k.txt
    $HaploC_sh -d ${wk_dir} -k $k >> $log_f 2>&1

    ############################## STEP4: call snps, do not depend on previous step

    k=snp
    log_f=$wk_dir/log_file/log.$k.txt
    $HaploC_sh -d ${wk_dir} -k $k >> $log_f 2>&1

    ############################## STEP5: intersect SNP and gzip

    k=sec
    log_f=$wk_dir/log_file/log.$k.txt
    for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk >> $log_f 2>&1 & done; wait

    ############################### STEP6: repair

    k=rep
    log_f=$wk_dir/log_file/log.$k.txt
    for chunk in $(seq 1 $n_chunks); do $HaploC_sh -d ${wk_dir} -k $k -n $chunk >> $log_f 2>&1 & done; wait

    ############################### STEP7: intg

    k=intg
    log_f=$wk_dir/log_file/log.$k.txt
    for chr in $(seq 1 23); do $HaploC_sh -d ${wk_dir} -k $k -c $chr >> $log_f 2>&1 & done; wait

    ############################### STEP8: opt

    k=opt
    log_f=$wk_dir/log_file/log.$k.txt
    $HaploC_sh -d ${wk_dir} -k $k >> $log_f 2>&1

    ############################### STEP9: mates

    k=mates
    log_f=$wk_dir/log_file/log.$k.txt
    $HaploC_sh -d ${wk_dir} -k $k >> $log_f 2>&1

    ###############################

    [ "$downstream" = "false" ] && exit

    k=diffIns
    log_f=$wk_dir/log_file/log.$k.txt
    bin_size=25000
    $downstream_sh -d $wk_dir -k $k -s $bin_size >> $log_f 2>&1

    k=diffComp
    log_f=$wk_dir/log_file/log.$k.txt
    $downstream_sh -d $wk_dir -k $k >> $log_f 2>&1

    k=HaploCNV
    log_f=$wk_dir/log_file/log.$k.txt
    bin_size=100000
    $downstream_sh -d $wk_dir -k $k -s $bin_size >> $log_f 2>&1
}

##############################

run_HaploC $ID &

##############################
