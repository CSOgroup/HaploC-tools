#!/bin/bash

#set -e ## This is causing problems; need better error detection
shopt -s extglob
export LC_ALL=C
juicer_version="1.13.02" 
### LOAD BWA AND SAMTOOLS

samtools=samtools
output=$($samtools 2>&1)

if [[ $output == *"error"* || $output == *"No such file"* ]]; then
  samtools=samtools ## conda intg
fi


SECONDS=0
######################################
add_mate_info() {

    local mate=$1
    if [[ "$mate" -eq 1 ]]
    then
        cat ${name1}.dedup.fastq | awk '{ if(NR%4==1) $1=$1 "_mate1"; print }' > $R12_with_mate
        rm  ${name1}.dedup.fastq
    elif [[ "$mate" -eq 2 ]]
    then
        cat ${name2}.dedup.fastq | awk '{ if(NR%4==1) $1=$1 "_mate2"; print }' > ${mate2_tmp}
        rm  ${name2}.dedup.fastq
    fi    
}

add_read_indicator() {
    local mate=$1
    if [[ "$mate" -eq 1 ]]
    then
        awk 'BEGIN{OFS="\t"}NF>=11{$1=$1"/1"; print}' ${name1}.n_sorted.sam > ${name1}.tmp.sam && mv ${name1}.tmp.sam ${name1}.n_sorted.sam
    elif [[ "$mate" -eq 2 ]]
    then
        awk 'BEGIN{OFS="\t"}NF>=11{$1=$1"/2"; print}' ${name2}.n_sorted.sam > ${name2}.tmp.sam && mv ${name2}.tmp.sam ${name2}.n_sorted.sam
    fi    
}

######################################

read1str="_R1" 
read2str="_R2" 

## Default options, overridden by command line arguments                
shortreadend=0
about=""
nofrag=1
justexact=0

while getopts "k:n:d:g:R:a:hrs:p:y:z:S:D:fjt:b:" opt; do
    case $opt in
    k) k=$OPTARG ;;        
    n) fastq_prefix=$OPTARG ;;        
    g) genomeID=$OPTARG ;;
    h) printHelpAndExit 0;;
    d) topDir=$OPTARG ;;
    s) site=$OPTARG ;;
    R) shortreadend=$OPTARG ;;
    r) shortread=1 ;;  #use short read aligner
    a) about=$OPTARG ;;
    p) genomePath=$OPTARG ;;  
    y) site_file=$OPTARG ;;
    z) refSeq=$OPTARG ;;
    S) stage=$OPTARG ;;
    D) juiceDir=$OPTARG ;;
    f) nofrag=0 ;; #use fragment maps
    b) ligation=$OPTARG ;;
        t) threads=$OPTARG ;;
    j) justexact=1 ;;
    [?]) printHelpAndExit 1;;
    esac
done


## If short read end is set, make sure it is 1 or 2
case $shortreadend in
    0) ;;
    1) ;;
    2) ;;
    *)  echo "$usageHelp"
    echo "$shortHelp2"
    exit 1
esac

threadstring="-t $threads"


## Directories to be created and regex strings for listing files
splitdir=${topDir}"/splits"
donesplitdir=$topDir"/done_splits"
fastqdir=${topDir}"/fastq"
outputdir=${topDir}"/aligned"

tmpdir=${topDir}"/HIC_tmp"

## https://stackoverflow.com/questions/8223170/bash-extracting-last-two-dirs-for-a-pathname/8223899
# grandparent=$(basename "$(dirname "$topDir")")/$(basename "$topDir")
# tmpdir=$TMPDIR/$grandparent/HIC_tmp ## this will not accelerate the time

logfile=${topDir}"/log.txt"
clumpy_logfile=${topDir}"/clumpy_log.txt"
rm2run_bwa=${topDir}"/rm2run_bwa.txt"
rm2run_pre=${topDir}"/rm2run_pre.txt"
bwa_thread_config=${topDir}"/replace_with_bwa_n.txt"


#################### set memory usage. I will use two cores

xGb=20G
if [[ $topDir == *"cso"* ]]; then
  xGb=20G
fi

echo -e "\nMemory requested: $xGb" >> $logfile

####################


## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]
then
    case $site in
    # HindIII) ligation="AAGCTAGCTT";;
    # DpnII) ligation="GATCGATC";;
    # MboI) ligation="GATCGATC";;
    # Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
    # Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
    # NcoI) ligation="CCATGCATGG";;

    HindIII) ligation="none";; ## I will not need to specify them here
    DpnII) ligation="none";;
    MboI) ligation="none";;
    Arima) ligation="none";;
    NcoI) ligation="none";;

    none) ligation="XXXX";;
    *)  ligation="XXXX"
        echo "$site not listed as recognized enzyme. Using $site_file as site file"
        echo "Ligation junction is undefined"
    esac
fi



echo $threads > ${bwa_thread_config}


## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

gzipped=1
read1='pseudo'
headfile=${outputdir}/header
name=$splitdir/${fastq_prefix}
fastq1=$fastqdir/${fastq_prefix}_R1.fastq.gz           
fastq2=$fastqdir/${fastq_prefix}_R2.fastq.gz           

name1=${name}_R1
name2=${name}_R2
usegzip=1

echo k:$k
# exit


if [[ "$k" == 0 ]]; then
{
    echo "I am in STEP_0"

    # rm -rf $outputdir
    # rm -rf $splitdir
    # rm -rf $tmpdir

    mkdir -p "$outputdir"
    mkdir -p "$splitdir"
    mkdir -p "$tmpdir"
    chmod 777 "$tmpdir"

    echo "rm2run_bwa" > ${rm2run_bwa}
    echo "rm2run_pre" > ${rm2run_pre}

    date > $headfile
    echo -ne 'Experiment description: ' >> $headfile
    echo -ne "Juicer version $juicer_version;" >> $headfile
    bwa 2>&1 | awk '$1=="Version:"{printf(" BWA %s; ", $2)}' >> $headfile 
    echo -ne "$threads threads; " >> $headfile
    java -version 2>&1 | awk 'NR==1{printf("%s; ", $0);}' >> $headfile 
    ${juiceDir}/scripts/juicer_tools -V 2>&1 | awk '$1=="Juicer" && $2=="Tools"{printf("%s; ", $0);}' >> $headfile
    echo "$0 $@" >> $headfile

    echo "Begin countligations at "`date` >> $logfile
    # source ${juiceDir}/scripts/common/countligations_pseudo.sh ## countligations.sh needs to gzip file and takes a lot of time. This value will not be write to the hic file / can check again
    echo -ne "7777777 " > ${name}_norm.txt.res.txt
    echo "7777777" > ${name}_linecount.txt

    ########################## remove duplicates
    echo "Begin clumpify at "`date` >> $logfile
    clumpify.sh -Xmx40g tmpdir=$fastqdir in=${fastq1} in2=${fastq2} out=${name1}.dedup.fastq.gz out2=${name2}.dedup.fastq.gz dedupe subs=5 t=20 2> ${clumpy_logfile}  ## https://www.biostars.org/p/225338/; default is subs=2, but I am more conserved
    echo "Finish clumpify at "`date` >> $logfile
    # exit

    echo "Finish STEP_0 at "`date` >> $logfile
    duration=$SECONDS
    echo "$(($duration / 3600)) hours or $(($duration / 60)) minutes elapsed."
    exit   
}
fi
     

if [[ "$k" == 1 ]]; then
{
    echo "I am in STEP_1"
    mkdir -p "$tmpdir"
    chmod 777 "$tmpdir"  
    rm ${fastq1} ${fastq2}
    touch ${fastq1} ${fastq2}

    ## ALIGN FASTQ AS SINGLE END, SORT BY READNAME, HANDLE CHIMERIC READS

    {
        ##########################
        echo "I am waiting you to remove at "`date` >> ${rm2run_bwa}
        while false; do
        if [[ ! -f "$rm2run_bwa" ]]; then
            break
        fi
        done

        ##########################
        threadstring=`cat ${bwa_thread_config}`

        bwa mem -t $threadstring $refSeq ${name1}.dedup.fastq.gz | sort --parallel=20 -S $xGb -T $tmpdir -k1,1f > ${name1}.n_sorted.sam ## potential issue if fastq name is ranked before SQ header
        # rm ${name1}.dedup.fastq.gz

        bwa mem -t $threadstring $refSeq ${name2}.dedup.fastq.gz | sort --parallel=20 -S $xGb -T $tmpdir -k1,1f > ${name2}.n_sorted.sam
        # rm ${name2}.dedup.fastq.gz

        $samtools view -@ 20 -bS ${name1}.n_sorted.sam -o ${name1}.n_sorted.bam
        $samtools view -@ 20 -bS ${name2}.n_sorted.sam -o ${name2}.n_sorted.bam

        echo "Finish bwa & sort at "`date` >> $logfile

        exit
    }
}
fi

if [[ "$k" == 2 ]]; then
{
    {
        echo "I am in STEP_2"
        mkdir -p "$tmpdir"
        chmod 777 "$tmpdir"  

        # samtools merge -@ 40 -n ${name}_R12.n_sorted.bam ${name1}.n_sorted.bam ${name2}.n_sorted.bam

        ########################## yliu, 11.04.2021
        # samtools sort ${name}_R12.n_sorted.bam -@ 40 -O BAM -o ${name}.pos_sorted.bam -m 2G
        # samtools index -@ 40 ${name}.pos_sorted.bam
        # bamCoverage -b ${name}.pos_sorted.bam --binSize 100 -p 40 --outFileName ${name}.pos_sorted.bw
        # echo "Finish bamCoverage at "`date` >> $logfile
        # # ##########################  
        
        # add read end indicator to readname
        for mate in 1 2; do add_read_indicator $mate & done
        wait 
        ##########################
        sort --parallel=20 -S $xGb -T $tmpdir -k1,1f -m ${name1}.n_sorted.sam ${name2}.n_sorted.sam > ${name}.n_sorted.sam
        rm ${name1}.n_sorted.sam ${name2}.n_sorted.sam
    }
    echo "Finish STEP_2 at "`date` >> $logfile

    duration=$SECONDS
    echo "$(($duration / 3600)) hours or $(($duration / 60)) minutes elapsed."
    exit
}
fi

if [[ "$k" == 3 ]]; then
{
    echo "I am in STEP_3"
    mkdir -p "$tmpdir"
    chmod 777 "$tmpdir"

    # call chimeric_blacklist.awk to deal with chimeric reads; 
    # sorted file is sorted by read name at this point
    touch ${name}_abnorm.sam ${name}_unmapped.sam
    awk -v "fname1"=${name}_norm.txt -v "fname2"=${name}_abnorm.sam -v "fname3"=${name}_unmapped.sam -f ${juiceDir}/scripts/common/chimeric_blacklist.awk ${name}.n_sorted.sam

    cut ${name}_norm.txt -f2,3,5,6 | awk '$1==$3' > ${name}.chr_pos
    awk ' $4 - $2 <= 10000 && $2 - $4 <= 10000 ' ${name}.chr_pos > ${name}.10kb.chr_pos
    awk ' $4 - $2 <= 5000 && $2 - $4 <= 5000 ' ${name}.10kb.chr_pos > ${name}.5kb.chr_pos
    awk ' $4 - $2 <= 1000 && $2 - $4 <= 1000 ' ${name}.10kb.chr_pos > ${name}.1kb.chr_pos
    mv ${name}.*kb.chr_pos $outputdir/

    exit

    rm $name.n_sorted.sam

    echo "Finish rm $name.sam "`date` >> $logfile
    printf "\n" >> $logfile

    # if any normal reads were written, find what fragment they correspond to 
    # and store that

    ${juiceDir}/scripts/common/fragment.pl ${name}_norm.txt $name.frag.txt $site_file                                                                
                                                                                    
    # sort by chromosome, fragment, strand, and position                    
    sort --parallel=20 -S $xGb -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $name.frag.txt > $name.sort.txt
    echo "Hi I do remove"
    rm ${name}_norm.txt ${name}.frag.txt

    #MERGE SORTED AND ALIGNED FILES
    sort --parallel=20 -S $xGb -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $splitdir/*.sort.txt  > $outputdir/merged_sort.txt
    rm $splitdir/*.sort.txt
    rm -r "$tmpdir"

    echo "Finish STEP_3 at "`date` >> $logfile
    duration=$SECONDS
    echo "$(($duration / 3600)) hours or $(($duration / 60)) minutes elapsed."
    exit
    
 }
fi

if [[ "$k" == 4 ]]; then
{
    #REMOVE DUPLICATES
    echo "I am in STEP_4"
    mkdir -p "$tmpdir"
    chmod 777 "$tmpdir"  

    touch ${outputdir}/dups.txt
    touch ${outputdir}/optdups.txt
    touch ${outputdir}/merged_nodups.txt
    
    if [ "$justexact" -eq 1 ]
    then
    awk -f ${juiceDir}/scripts/common/dups.awk -v name=${outputdir}/ -v nowobble=1 ${outputdir}/merged_sort.txt
    else
    awk -f ${juiceDir}/scripts/common/dups.awk -v name=${outputdir}/ ${outputdir}/merged_sort.txt
    fi

    mv ${outputdir}/optdups.txt ${outputdir}/opt_dups.txt 
    rm ${outputdir}/merged_sort.txt

    export _JAVA_OPTIONS=-Xmx16384m
    export LC_ALL=en_US.UTF-8 
    tail -n1 $headfile | awk '{printf"%-1000s\n", $0}' > $outputdir/inter.txt;
    cat $splitdir/*.res.txt | awk -f ${juiceDir}/scripts/common/stats_sub.awk >> $outputdir/inter.txt
    ${juiceDir}/scripts/common/juicer_tools LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt
    # cp $outputdir/inter.txt $outputdir/inter_30.txt 

    # ${juiceDir}/scripts/common/statistics.pl -s $site_file -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt 
    ${juiceDir}/scripts/common/statistics.pl -s "none" -l "none" -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt 

    cat $splitdir/*_abnorm.sam > $outputdir/abnormal.sam
    cat $splitdir/*_unmapped.sam > $outputdir/unmapped.sam
    awk -f ${juiceDir}/scripts/common/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
    # Collisions dedupping script goes here
    echo "Finish STEP_4 at "`date` >> $logfile

    rm $outputdir/abnormal.sam $outputdir/unmapped.sam $outputdir/dups.txt $outputdir/opt_dups.txt $outputdir/collisions.txt
    rm ${name}_abnorm.sam ${name}_unmapped.sam
    rm ${name}_norm.txt.res.txt
    rm ${name}_linecount.txt

    duration=$SECONDS
    echo "$(($duration / 3600)) hours or $(($duration / 60)) minutes elapsed."
    exit
}
fi

if [[ "$k" == 5 ]]; then
{
    #STATISTICS
    #Skip if post-processing only is required
        
    echo "I am waiting you to remove at "`date` >> ${rm2run_pre}
    while false; do
    if [[ ! -f "$rm2run_pre" ]]; then
        break
    fi
    done

    #CREATE HIC FILES
    # if early exit, we stop here, once the statistics are calculated

    ${juiceDir}/scripts/common/juicer_tools pre -r 10000,25000,40000,50000,100000,500000,1000000 -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/inter.hic $genomePath

    #CHECK THAT PIPELINE WAS SUCCESSFUL
    export early=$earlyexit
    export splitdir=$splitdir
}
fi
