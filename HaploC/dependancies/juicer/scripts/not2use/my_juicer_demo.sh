

juicer_sh=/mnt/ndata/Yuanlong/3.Packages/juicer/scripts/my_juicer.sh
k=0 ## I split the whole thing into 6 steps // k=0: use clumpify to depulicate the fastq files; k=1: bwa alignment; k=2: adding mate infomation into the bam/sam file; k=3: generate merged_sort.txt; k=4: generate merged_nodups.txt; k=5: create hic file
enzyme_name=MboI
wk_dir=/mnt/ptemp/Yuanlong/unil_nas/default_nonsensitive/D2c/Yuanlong//1.Projects/15.Phasing_SV/1.Data/HCT_116/run10 ## the working directory for this hic study, should contain the fastq file under ./fastq
fastq_prefix=SRR6107791 ## the prefix of fastq, should be named as ${fastq_prefix}_R1.fastq.gz and ${fastq_prefix}_R2.fastq.gz
threads=10
log_file=xx


$juicer_sh -k $k -s $enzyme_name -n $fastq_prefix -g hg19 -z /mnt/ndata/Yuanlong/3.Packages/juicer/references/hg19_chr1_22XYM.fa -d $wk_dir -D /mnt/ndata/Yuanlong/3.Packages/juicer -y /mnt/ndata/Yuanlong/3.Packages/juicer/restriction_sites/hg19_${enzyme_name}.txt -p /mnt/ndata/Yuanlong/3.Packages/juicer/references/hg19.chrom.sizes -t $threads >> $log_file