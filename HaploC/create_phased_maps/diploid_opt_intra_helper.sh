	#!/bin/bash
	
	wk_dir=$1
	impute=$2 ## impute or n_impute
	repo_dir=$3
	genome_version=$4

	##########################################################

	identity=opt
	wk_dir_mega=$wk_dir/mega/
	merged_nodups_snp_dir=$wk_dir_mega/aligned/merged_nodup_chrXchrY

	mkdir -p $wk_dir_mega/HapCut/hic_phased/

 	##########################################################
   
    juicerDir=${repo_dir}/HaploC/dependancies/juicer/
    juicer_jar=${juicerDir}/scripts/juicer_tools_1.19.02.jar
	htrans_cis_info_global_fun=${repo_dir}/HaploC/H-trans-evaluation/htrans_ratio_global.R
	# htrans_cis_info_local_fun=${repo_dir}/HaploC/htrans_ratio_local.R

    #########################################################

    chr2keep_f=$wk_dir/mega/VCFs/chr2phase.txt
	chr_nums=($(<$chr2keep_f))

	#########################################################

	helper_for_chr()
	{
		local chr_num_x=$1
		echo $chr_num_x

		chr_num_y=$chr_num_x
		chr_num=$chr_num_x

		merged_nodups_f=${merged_nodups_snp_dir}/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}.snps

		chr_pos_f=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num}_$identity.chr_pos
		snp_pat_mat_f=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num}_$identity.pat_mat

		if [ $impute == "n_impute" ]
		then
		   	prefix=$wk_dir_mega/HapCut/hic_phased/hic_chr${chr_num}_$identity.n_impute
		else
			prefix=$wk_dir_mega/HapCut/hic_phased/hic_chr${chr_num}_$identity
		fi

		dip_out_f=${prefix}.dip_out

		
		rm -f ${dip_out_f}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.mat.txt
		rm -f ${prefix}.pat.txt
		rm -f ${prefix}.diff.txt

		echo $snp_pat_mat_f

		awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' $merged_nodups_f | $juicerDir/scripts/common/diploid.pl -s $chr_pos_f -o $snp_pat_mat_f > $dip_out_f	
		awk -v stem=$prefix. -f ${juicerDir}/scripts/common/diploid_split.awk $dip_out_f

		# exit

		######################################################### get htrans profile, both locally and globally

	    slim_out=${dip_out_f}.slim
	   
	    cut -f3,7,12,13 $dip_out_f -d" " | grep -v "mismatch" | awk 'length($3)>4 && length($4)>4' | cut -f3,4 -d" " > $slim_out
	    wc -l $slim_out

	  	if [[ $impute == "impute" ]]
		then
			Rscript $htrans_cis_info_global_fun $wk_dir $identity $chr_num_x generate_data
			# Rscript $htrans_cis_info_local_fun $wk_dir_mega $identity $chr_num_x
		fi  


		# rm $slim_out
		# awk -v stem=$prefix. -f ${juicerDir}/scripts/common/diploid_split_inter_YL.awk $dip_out_f

		## with imputation: 

		if [[ $impute == "n_impute" ]]
		then
		   	sort -k2,2d -m ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.mat.txt 
			sort -k2,2d -m ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.pat.txt
		else
		   	sort -k2,2d -m ${prefix}.maternal.txt ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.mat.txt 
			sort -k2,2d -m ${prefix}.paternal.txt ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.pat.txt			
		fi

		echo ${prefix}.pat.txt

		rm -f ${dip_out_f}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.diff.txt
	}


	#########################################################

	generate_hic_impute()
	{
		local wk_dir_mega=$1
		local which=$2
		local genome_version=$3

		# cat $wk_dir_mega/HapCut/hic_phased/hic_chr*$identity.$which.txt > $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.txt

		to_combine=$(printf " $wk_dir_mega/HapCut/hic_phased/hic_chr%s_$identity.$which.txt" "${chr_nums[@]}")
		cat $to_combine > $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.txt
	
		if [ -z "$hic_name" ]
		then
		    java -Xmx200g -jar $juicer_jar pre -r 10000,20000,25000,40000,50000,60000,70000,80000,100000,250000,500000,1000000 $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.txt $wk_dir_mega/HapCut/hic_phased/all_chrs.$identity.$which.hic $genome_version
		else
		   	java -Xmx200g -jar $juicer_jar pre -r 10000,20000,25000,40000,50000,60000,70000,80000,100000,250000,500000,1000000 $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.txt $wk_dir_mega/HapCut/hic_phased/all_chrs.$identity.$hic_name.$which.hic $genome_version
		fi

	}

	generate_hic_n_impute()
	{
		local wk_dir_mega=$1
		local which=$2
		local genome_version=$3

		# cat $wk_dir_mega/HapCut/hic_phased/hic_chr*$identity.$which.txt > $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.txt

		to_combine=$(printf " $wk_dir_mega/HapCut/hic_phased/hic_chr%s_$identity.n_impute.$which.txt" "${chr_nums[@]}")
		cat $to_combine > $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.n_impute.$which.txt
	
		if [ -z "$hic_name" ]
		then
		    java -Xmx200g -jar $juicer_jar pre -r 10000,20000,25000,40000,50000,60000,70000,80000,100000,250000,500000,1000000 $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.n_impute.$which.txt $wk_dir_mega/HapCut/hic_phased/all_chrs.$identity.n_impute.$which.hic $genome_version
		else
		   	java -Xmx200g -jar $juicer_jar pre -r 10000,20000,25000,40000,50000,60000,70000,80000,100000,250000,500000,1000000 $wk_dir_mega/HapCut/hic_phased/hic_ALL_chrs.$identity.n_impute.$which.txt $wk_dir_mega/HapCut/hic_phased/all_chrs.$identity.n_impute.$hic_name.$which.hic $genome_version
		fi

	}


	#########################################################

	for chr_num_x in ${chr_nums[@]} ; do helper_for_chr $chr_num_x & done
	# for chr_num_x in ${chr_nums[@]} ; do helper_for_chr $chr_num_x; done	
	wait

	#########################################################

	if [ $impute == "n_impute" ]
	then
		for which in mat pat; do generate_hic_n_impute $wk_dir_mega $which $genome_version & done
	else
		for which in mat pat; do generate_hic_impute $wk_dir_mega $which $genome_version & done
		Rscript $htrans_cis_info_global_fun $wk_dir $identity psudo_chr generate_plot		
	fi
	