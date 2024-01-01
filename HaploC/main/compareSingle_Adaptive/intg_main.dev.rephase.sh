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

		############################## STEP2: align and so on	

		k=bwa
		log_f=$wk_dir/log_file/log.$k
		n_chunks=$(cat "$n_chunks_f")
		# jid_bwa=$(sbatch -J $k.$line -c $thread4bwa --mem 120G --time=300 --array=1-$n_chunks -o $log_f.chunk=%a.txt $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## for files > 10Gb
		echo bwa $jid_bwa >> $job_id_f
		# exit

		############################## STEP3: proceed to generate hic

		k=hic
		log_f=$wk_dir/log_file/log.$k
		# jid_hic=$(sbatch -J $k.$line -c 1 --mem 80G --time=300 --array=1-$n_chunks -o $log_f.chunk=%a.txt --dependency=afterany:$jid_bwa $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		echo hic $jid_hic >> $job_id_f
		# exit

		############################## STEP3x: generate mega / do not do mega

		k=mega
		log_f=$wk_dir/log_file/log.$k.txt
		# jid_mega=$(sbatch -J $k.$line -c 1 --mem 100G --time=1000 -o $log_f --dependency=afterany:$jid_hic $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		echo mega $jid_mega >> $job_id_f
		# exit

		############################## STEP4: call snps, do not depend on previous step

		k=snp 
		log_f=$wk_dir/log_file/log.$k.txt
		# jid_snp=$(sbatch -J $k.$line -c 20 --mem 100G --time=1000 -o $log_f --dependency=afterany:$jid_bwa $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
		echo snp $jid_snp >> $job_id_f

		############################## STEP5: intersect SNP and gzip

		k=sec 
		log_f=$wk_dir/log_file/log.$k
		# jid_sec=$(sbatch -J $k.$line -c 1 --mem 40G --time=300 --array=1-$n_chunks -o $log_f.chunk=%a.txt --dependency=afterany:$jid_hic:$jid_snp $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		echo sec $jid_sec >> $job_id_f	

		# exit	

		############################### STEP6: repair

		k=rep
		log_f=$wk_dir/log_file/log.$k
		# jid_rep=$(sbatch -J $k.$line -c 23 --mem 100G --time=600 --array=1-$n_chunks -o $log_f.chunk=%a.txt --dependency=afterany:$jid_bwa $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
		echo rep $jid_rep >> $job_id_f		

		############################### STEP7: intg

		k=intg
		log_f=$wk_dir/log_file/log.$k
		# jid_intg=$(sbatch -J $k.$line -c 7 --mem 70G --time=1000 --array=1-23 -o $log_f.chr=%a.txt $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		echo intg $jid_intg >> $job_id_f	

		# exit	

		############################### STEP8: opt

		k=opt
		log_f=$wk_dir/log_file/log.$k.txt
		jid_opt=$(sbatch -J $k.$line -c 23 --mem 100G --time=100 -o $log_f $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		echo opt $jid_opt >> $job_id_f		

		############################### STEP9: mates

		k=mates
		log_f=$wk_dir/log_file/log.$k.txt
		jid_mates=$(sbatch -J $k.$line -c 23 --mem 150G --time=600 -o $log_f --dependency=afterany:$jid_opt $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		echo mates $jid_mates >> $job_id_f

		k=htrans
		log_f=$wk_dir/log_file/log.$k.txt
		jid_htrans=$(sbatch -J $k.$line -c 23 --mem 100G --time=100 -o $log_f --dependency=afterany:$jid_mates $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
		echo htrans $jid_htrans >> $job_id_f

		k=htransplot
		log_f=$wk_dir/log_file/log.$k.txt
		jid_htransplot=$(sbatch -J $k.$line -c 1 --mem 50G --time=100 -o $log_f --dependency=afterany:$jid_htrans $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
		echo htransplot $jid_htransplot >> $job_id_f

	}

	##############################

	HaploC $line &

	##############################