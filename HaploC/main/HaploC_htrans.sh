	#!/bin/bash
	export LC_ALL=C

	##############################

	# find ./ -type d -name "*REDO" -exec env unil_script_dir="$unil_script_dir" sh -c '$unil_script_dir/HaploC.htrans.sh "$(basename "{}")"' \;

	unil_work=/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/
	script_dir=$unil_work/0.Scripts/HaploC/main/
	phasing_dir=/scratch/yliu5/1.Projects/15.Phasing/1.Data/

	##############################

	line=$1
	wk_dir=$phasing_dir/$line
	echo $line

	##############################

	job_id_f=$wk_dir/log_file/job_ids.txt

	############################## STEP0: split fastq

	job_sh=$script_dir/intg_jobs.ALL.sh

	##############################

	HaploC()
	{
		unil_work=/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/
		script_dir=$unil_work/0.Scripts/HaploC/main/
		phasing_dir=/scratch/yliu5/1.Projects/15.Phasing/1.Data/
		job_sh=$script_dir/intg_jobs.ALL.sh

		line=$1
		wk_dir=$phasing_dir/$line

		############################### STEP10: htrans

		# k=opt
		# log_f=$wk_dir/log_file/log.$k.txt
		# jid_opt=$(sbatch -J $k.$line -c 23 --mem 100G --time=150 -o $log_f $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		# echo opt $jid_opt >> $job_id_f

		# ############################### STEP9: mates

		k=htrans
		log_f=$wk_dir/log_file/log.$k.txt
		jid_htrans=$(sbatch -J $k.$line -c 23 --mem 100G --time=100 -o $log_f $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
		echo htrans $jid_htrans >> $job_id_f

		k=htransplot
		log_f=$wk_dir/log_file/log.$k.txt
		jid_htransplot=$(sbatch -J $k.$line -c 1 --mem 50G --time=100 -o $log_f --dependency=afterany:$jid_htrans $job_sh $wk_dir $k | sed 's/Submitted batch job //') ## memory setting should be sufficient as in the free script / samtools merge
		echo htransplot $jid_htransplot >> $job_id_f

		k=mates
		log_f=$wk_dir/log_file/log.$k.txt
		jid_mates=$(sbatch -J $k.$line -c 23 --mem 250G --time=600 -o $log_f --dependency=afterany:$jid_htransplot $job_sh $wk_dir $k | sed 's/Submitted batch job //')
		echo mates $jid_mates >> $job_id_f

	}

	##############################

	HaploC $line &

	##############################