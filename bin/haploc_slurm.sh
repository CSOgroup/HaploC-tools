#!/bin/bash --login
export LC_ALL=C

##############################

usage() {
    echo "Usage: haploc_slurm.sh [-d wk_dir] [-x] [-h]"
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

config_f=$wk_dir/config/config.txt
ID=$(basename "$wk_dir"); echo $ID
sizeGb=$(grep '^sizeGb=' $config_f | cut -d '=' -f 2)
thread4bwa=$(grep '^thread4bwa=' $config_f | cut -d '=' -f 2)
repo_dir=$(grep '^repo_dir=' $config_f | cut -d '=' -f 2)
HaploC_sh=$repo_dir/HaploC/main/HaploC_helper_slurm.sh
downstream_sh=$repo_dir/bin/downstreams.sh

##############################

mkdir -p $wk_dir/log_file
echo $wk_dir
job_id_f=$wk_dir/log_file/job_ids.txt

run_HaploC()
{
    ############################## STEP0: split fastq
    
    k=split
    log_f=$wk_dir/log_file/log.$k
    data_source_f=$wk_dir/log_file/data_source.txt
    n_data=$(find $wk_dir -type f -name '*1.fastq.gz' | wc -l)
    find $wk_dir -type f -name '*_1.fastq.gz' -exec basename {} \; | sed 's/_1.fastq.gz//' > $data_source_f
    jid_split=$(sbatch -J split.$ID -c 5 --mem 50G --time=100 --array=1-$n_data -o $log_f.data=%a.txt $HaploC_sh -d $wk_dir -k $k $sizeGb | sed 's/Submitted batch job //') ## 5 is set as the number of threads for seqkit. More threads seems not increase the performance
    echo split $jid_split > $job_id_f
    
    ############################## STEP1: distribute data to runs
    
    k=distr
    log_f=$wk_dir/log_file/log.$k.txt
    jid_distr=$(sbatch -J distr.$ID -c 1 --mem 10G --time=30 -o $log_f --dependency=afterany:$jid_split $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //')
    
    ##############################
    
    n_chunks_f=$wk_dir/log_file/n_chunks.txt
    while [ ! -f $n_chunks_f ]; do sleep 10; done
    
    ############################## STEP2: align and so on
    
    k=bwa;
    log_f=$wk_dir/log_file/log.$k
    n_chunks=$(cat "$n_chunks_f")
    jid_bwa=$(sbatch -J $k.$ID -c $thread4bwa --mem 120G --time=300 --array=1-$n_chunks -o $log_f.chunk=%a.txt $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //') ## for files > 10Gb
    echo bwa $jid_bwa >> $job_id_f
    # exit
    
    ############################## STEP3: proceed to generate hic
    
    k=hic
    log_f=$wk_dir/log_file/log.$k
    jid_hic=$(sbatch -J $k.$ID -c 1 --mem 80G --time=300 --array=1-$n_chunks -o $log_f.chunk=%a.txt --dependency=afterany:$jid_bwa $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //')
    echo hic $jid_hic >> $job_id_f
    # exit
    
    ############################## STEP3x: generate mega / do not do mega
    
    k=mega
    log_f=$wk_dir/log_file/log.$k.txt
    jid_mega=$(sbatch -J $k.$ID -c 1 --mem 100G --time=1000 -o $log_f --dependency=afterany:$jid_hic $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //')
    echo mega $jid_mega >> $job_id_f
    
    ############################## STEP4: call snps, do not depend on previous step
    
    k=snp
    log_f=$wk_dir/log_file/log.$k.txt
    jid_snp=$(sbatch -J $k.$ID -c 20 --mem 100G --time=1000 -o $log_f --dependency=afterany:$jid_bwa $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
    echo snp $jid_snp >> $job_id_f
    
    ############################## STEP5: intersect SNP and gzip
    
    k=sec
    log_f=$wk_dir/log_file/log.$k
    jid_sec=$(sbatch -J $k.$ID -c 1 --mem 20G --time=300 --array=1-$n_chunks -o $log_f.chunk=%a.txt --dependency=afterany:$jid_snp:$jid_hic $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //')
    echo sec $jid_sec >> $job_id_f
    
    ############################### STEP6: repair
    
    k=rep
    log_f=$wk_dir/log_file/log.$k
    jid_rep=$(sbatch -J $k.$ID -c 5 --mem 50G --time=100 --array=1-$n_chunks -o $log_f.chunk=%a.txt --dependency=afterany:$jid_bwa $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
    echo rep $jid_rep >> $job_id_f
    
    ############################### STEP7: intg
    
    k=intg
    log_f=$wk_dir/log_file/log.$k
    jid_intg=$(sbatch -J $k.$ID -c 7 --mem 50G --time=1000 --array=1-23 -o $log_f.chr=%a.txt --dependency=afterany:$jid_rep:$jid_snp $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //')
    echo intg $jid_intg >> $job_id_f
    
    ############################### STEP8: opt
    
    k=opt
    log_f=$wk_dir/log_file/log.$k.txt
    jid_opt=$(sbatch -J $k.$ID -c 23 --mem 100G --time=100 -o $log_f --dependency=afterany:$jid_intg:$jid_sec $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //')
    echo opt $jid_opt >> $job_id_f
    
    ############################### STEP9: mates
    
    k=mates
    log_f=$wk_dir/log_file/log.$k.txt
    jid_mates=$(sbatch -J $k.$ID -c 23 --mem 50G --time=600 -o $log_f --dependency=afterany:$jid_opt $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //')
    echo mates $jid_mates >> $job_id_f
    
    # k=htrans
    # log_f=$wk_dir/log_file/log.$k.txt
    # jid_htrans=$(sbatch -J $k.$ID -c 23 --mem 100G --time=100 -o $log_f --dependency=afterany:$jid_mates $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
    # echo htrans $jid_htrans >> $job_id_f
    
    # k=htransplot
    # log_f=$wk_dir/log_file/log.$k.txt
    # jid_htransplot=$(sbatch -J $k.$ID -c 1 --mem 50G --time=100 -o $log_f --dependency=afterany:$jid_htrans $HaploC_sh -d $wk_dir -k $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
    # echo htransplot $jid_htransplot >> $job_id_f
    
    ##############################
    
    [ "$downstream" = "false" ] && exit
    
    k=diffIns
    bin_size=25000
    log_f=$wk_dir/log_file/log.$k.txt
    jid_diffIns=$(sbatch -J $k.$ID -c 10 --mem 50G --time=30 -o $log_f --dependency=afterany:$jid_mates $downstream_sh -d $wk_dir -k $k -s $bin_size | sed 's/Submitted batch job //') ## for files > 10Gb
    echo diffIns $jid_diffIns >> $job_id_f
    
    k=diffComp
    log_f=$wk_dir/log_file/log.$k.txt
    jid_diffComp=$(sbatch -J $k.$ID -c 10 --mem 50G --time=30 -o $log_f --dependency=afterany:$jid_mates $downstream_sh -d $wk_dir -k $k | sed 's/Submitted batch job //') ## for files > 10Gb
    echo diffComp $jid_diffComp >> $job_id_f
    
    k=HaploCNV
    bin_size=100000
    log_f=$wk_dir/log_file/log.$k.txt
    jid_HaploCNV=$(sbatch -J $k.$ID -c 10 --mem 50G --time=30 -o $log_f --dependency=afterany:$jid_diffComp $downstream_sh -d $wk_dir -k $k -s $bin_size | sed 's/Submitted batch job //')
    echo HaploCNV $jid_HaploCNV >> $job_id_f
}

##############################

run_HaploC $ID &

##############################