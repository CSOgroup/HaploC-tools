	#!/bin/bash
	export LC_ALL=C

	##############################

	unil_work=/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/
	script_dir=$unil_work/0.Scripts/HaploC/main/
	phasing_dir=/scratch/yliu5/1.Projects/15.Phasing/1.Data/

	##############################

	line=$1
	genome_version=hg19
	wk_dir=$phasing_dir/$line
	echo $line

	sizeGB=8
	thread4bwa=30

	##############################

	mkdir -p $wk_dir/log_file
	echo $wk_dir
	job_id_f=$wk_dir/log_file/job_ids.txt

	##############################

	STORAGE_USE()
	{
		unil_work=/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/
		script_dir=$unil_work/0.Scripts/HaploC/main/
		phasing_dir=/scratch/yliu5/1.Projects/15.Phasing/1.Data/

		line=$1
		wk_dir=$phasing_dir/$line
		Rscript $script_dir/storage_use.R $wk_dir
	}

	# STORAGE_USE $line &

	############################## STEP0: split fastq

	k=split
	log_f=$wk_dir/log_file/log.$k
	job_sh=$script_dir/intg_jobs.ALL.sh
	SRA_data_f=$unil_work/0.Scripts/SRA_database/SRA.$line.tsv
	n_data=$(wc -l < "$SRA_data_f") ## number of fastqs
	# jid_split=$(sbatch -J split.$line -c 5 --mem 100G --time=500 --array=1-$n_data -o $log_f.data=%a.txt $job_sh $wk_dir $k $sizeGB | sed 's/Submitted batch job //') ## 5 is set as the number of threads for seqkit. More threads seems not increase the performance
	echo split $jid_split > $job_id_f
		
	############################## STEP1: distribute data to runs

	k=org
	job_sh=$script_dir/intg_jobs.ALL.sh
	# jid_org=$(sbatch -J org.$line -c 1 --mem 10G --time=30 --dependency=afterany:$jid_split $job_sh $wk_dir $k | sed 's/Submitted batch job //')

	##############################

	HaploC()
	{
		unil_work=/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/
		script_dir=$unil_work/0.Scripts/HaploC/main/
		phasing_dir=/scratch/yliu5/1.Projects/15.Phasing/1.Data/
		job_sh=$script_dir/intg_jobs.ALL.sh

		line=$1
		wk_dir=$phasing_dir/$line
		echo HaploC for $wk_dir
		n_chunks_f=$wk_dir/log_file/n_chunks.txt

		##############################

		echo "WAITING HaploC for $line" > $wk_dir/step2to5_status.txt
		while [ ! -f $n_chunks_f ]
		do
		  sleep 10
		done
		echo "HaploC STARTED" >> $wk_dir/step2to5_status.txt

		############################## STEP4: call snps, do not depend on previous step

		k=snp 
		log_f=$wk_dir/log_file/log.$k.txt
		jid_snp=$(sbatch -J $k.$line -c 20 --mem 100G --time=1000 -o $log_f $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
		echo snp $jid_snp >> $job_id_f
	}

	##############################

	HaploC $line &

	##############################