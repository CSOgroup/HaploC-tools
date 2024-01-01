	#!/bin/bash

	## construct diploid hic map, for inter
	## yliu, 16-04-2021
	wk_dir_mega=$1
	chr_num_x=$2
	chr_num_y=$3
	identity=$4
	create_hic=$5

	impute=$6 ## impute or n_impute
	repo_dir=$7
	genome_version=$8

	juicerDir=${repo_dir}/HaploC/dependancies/juicer/
	juicer_jar=${juicerDir}/scripts/juicer_tools_1.19.02.jar
		
	mkdir -p $wk_dir_mega/HapCut/hic_phased/

	######################################## 

	if [[ "$create_hic" == "no" ]]; then

		merged_nodups_f=$wk_dir_mega/aligned/merged_nodup_chrXchrY/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}.snps

		chr_pos_f_x=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num_x}_$identity.chr_pos
		chr_pos_f_y=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num_y}_$identity.chr_pos

		snp_pat_mat_f_x=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num_x}_$identity.pat_mat
		snp_pat_mat_f_y=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num_y}_$identity.pat_mat
		
		chr_pos_f_xy=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num_x}chr${chr_num_y}_$identity.chr_pos
		snp_pat_mat_f_xy=$wk_dir_mega/HapCut/phased_hap_intg/hic_chr${chr_num_x}chr${chr_num_y}_$identity.pat_mat

		cat $chr_pos_f_x $chr_pos_f_y > $chr_pos_f_xy
		cat $snp_pat_mat_f_x $snp_pat_mat_f_y > $snp_pat_mat_f_xy

		if [[ $impute == "n_impute" ]]
		then
		   	prefix=$wk_dir_mega/HapCut/hic_phased/hic_chr${chr_num_x}chr${chr_num_y}_$identity.n_impute
		   	prefix_slim=$wk_dir_mega/HapCut/hic_phased/hic_interALL_$identity.n_impute
		else
			prefix=$wk_dir_mega/HapCut/hic_phased/hic_chr${chr_num_x}chr${chr_num_y}_$identity
		   	prefix_slim=$wk_dir_mega/HapCut/hic_phased/hic_interALL_$identity			
		fi

		dip_out_f=${prefix}.dip_out


		# awk -v chr=$i '$2==chr && $6==chr && $9 >= 10 && $12 >= 10' $merged_nodups_f | ${juiceDir}/scripts/diploid.pl -s ${stem}${i}_chr_pos.txt -o ${stem}${i}_paternal_maternal.txt > diploid_${i}.txt
		rm -f ${dip_out_f}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.mat.txt
		rm -f ${prefix}.pat.txt
		rm -f ${prefix}.diff.txt

		awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' $merged_nodups_f | $juicerDir/scripts/common/diploid.pl -s $chr_pos_f_xy -o $snp_pat_mat_f_xy > $dip_out_f
		awk -v stem=$prefix. -f ${juicerDir}/scripts/common/diploid_split.awk $dip_out_f

		# sort -k2,2d -m ${prefix}.maternal.txt ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.mat.txt 
		# sort -k2,2d -m ${prefix}.paternal.txt ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.pat.txt 

		# ## only for the real translocated chr-pair, imputation can be used

		if [[ $impute == "n_impute" ]]
		then
			sort -k2,2d -m ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.mat.txt 
			sort -k2,2d -m ${prefix}.paternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.pat.txt 
		else
			sort -k2,2d -m ${prefix}.maternal.txt ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.mat.txt 
			sort -k2,2d -m ${prefix}.paternal.txt ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.pat.txt	
		fi		

		rm -f ${dip_out_f}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.mat.txt
		rm -f ${prefix}.pat.txt
		rm -f ${prefix}.diff.txt
		rm -f ${prefix}.mismatch.txt

		########################################### switch mat / pat of chr_num_y: switch.pat: pat vs mat | switch.mat: mat vs pat

		chr_pos_f_xy_switch=$chr_pos_f_xy
		snp_pat_mat_f_xy_switch=${snp_pat_mat_f_xy}_switch
		awk '{print($1,$2,$3)}' ${snp_pat_mat_f_x} > $snp_pat_mat_f_xy_switch
		awk '{print($1,$3,$2)}' ${snp_pat_mat_f_y} >> $snp_pat_mat_f_xy_switch


		if [[ $impute == "n_impute" ]]
		then
		   	prefix=$wk_dir_mega/HapCut/hic_phased/hic_chr${chr_num_x}chr${chr_num_y}_$identity.n_impute.switch
		   	prefix_slim=$wk_dir_mega/HapCut/hic_phased/hic_interALL_$identity.n_impute.switch
		else
			prefix=$wk_dir_mega/HapCut/hic_phased/hic_chr${chr_num_x}chr${chr_num_y}_$identity.switch
		   	prefix_slim=$wk_dir_mega/HapCut/hic_phased/hic_interALL_$identity.switch			
		fi

		dip_out_f=${prefix}.dip_out

		###########################################

		rm -f ${dip_out_f}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.mat.txt
		rm -f ${prefix}.pat.txt
		rm -f ${prefix}.diff.txt
		rm -f ${prefix}.mismatch.txt

		# awk -v chr=$i '$2==chr && $6==chr && $9 >= 10 && $12 >= 10' $merged_nodups_f | ${juiceDir}/scripts/diploid.pl -s ${stem}${i}_chr_pos.txt -o ${stem}${i}_paternal_maternal.txt > diploid_${i}.txt
		
		awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' $merged_nodups_f | $juicerDir/scripts/common/diploid.pl -s $chr_pos_f_xy_switch -o $snp_pat_mat_f_xy_switch > $dip_out_f
		awk -v stem=$prefix. -f ${juicerDir}/scripts/common/diploid_split.awk $dip_out_f
		
		## only for the real translocated chr-pair, imputation can be used

		if [[ $impute == "n_impute" ]]
		then
			sort -k2,2d -m ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.mat.txt 
			sort -k2,2d -m ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.pat.txt 
		else
			sort -k2,2d -m ${prefix}.maternal.txt ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.mat.txt 
			sort -k2,2d -m ${prefix}.paternal.txt ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' >> ${prefix_slim}.pat.txt 
		fi


		rm -f ${dip_out_f}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.mat.txt
		rm -f ${prefix}.pat.txt
		rm -f ${prefix}.diff.txt
		rm -f ${prefix}.mismatch.txt

		rm $snp_pat_mat_f_xy
		rm $chr_pos_f_xy_switch $snp_pat_mat_f_xy_switch
		echo ${prefix_slim}.pat.txt 
	fi

	###########################################

	generate_hic_impute()
	{
		local wk_dir_mega=$1
		local which=$2
		local genome_version=$3
		text_f=$wk_dir_mega/HapCut/hic_phased/hic_interALL_$identity.$which.txt
		echo $text_f
		java -Xmx200g -jar $juicer_jar pre -r 10000,20000,25000,40000,50000,100000,500000,1000000 $text_f $wk_dir_mega/HapCut/hic_phased/all_chrs.$identity.inter.$which.hic $genome_version
	}

	###########################################

	generate_hic_n_impute()
	{
		local wk_dir_mega=$1
		local which=$2
		local genome_version=$3		
		text_f=$wk_dir_mega/HapCut/hic_phased/hic_interALL_$identity.n_impute.$which.txt
		echo $text_f
		java -Xmx200g -jar $juicer_jar pre -r 10000,20000,25000,40000,50000,100000,500000,1000000 $text_f $wk_dir_mega/HapCut/hic_phased/all_chrs.$identity.inter.n_impute.$which.hic $genome_version
	}

	###########################################

	if [[ "$create_hic" == "yes" ]]; then
		echo "I am creating hap inter hic files"
		if [[ $impute == "n_impute" ]]
		then
			for which in mat pat switch.mat switch.pat; do generate_hic_n_impute $wk_dir_mega $which $genome_version & done
		else
			for which in mat pat switch.mat switch.pat; do generate_hic_impute $wk_dir_mega $which $genome_version & done
		fi
	fi

	###########################################