	#!/bin/bash

	## construct diploid hic map, for inter
	wk_dir=$1
	repo_dir=$2
	genome_version=$3

 	##########################################################
 	
 	wk_dir_mega=$wk_dir/mega
    juicerDir=${repo_dir}/HaploC/dependancies/juicer/
    juicer_jar=${juicerDir}/scripts/juicer_tools_1.19.02.jar

	########################################

	hic_phased_dir=$wk_dir_mega/HapCut/hic_phased/
	mkdir -p $hic_phased_dir/generated_random_mix
	cd $hic_phased_dir/generated_random_mix

	awk '{print > "hic_"$3".opt.mat.txt" }' $hic_phased_dir/hic_ALL_chrs.opt.mat.txt 
	awk '{print > "hic_"$3".opt.pat.txt" }' $hic_phased_dir/hic_ALL_chrs.opt.pat.txt 


	generate_hic()
	{
		local which=$1
		local k=$2
		local genome_version=$3
		text=$hic_phased_dir/generated_random_mix/all_chrs.opt.k=$k.tmp$which.txt
		hic=$hic_phased_dir/generated_random_mix/all_chrs.opt.k=$k.rand$which.hic
		java -Xmx200g -jar $juicer_jar pre -r 20000,40000,50000,100000,250000,500000,1000000 $text $hic $genome_version
	}

	#########################################################

	chr2keep_f=$wk_dir/mega/VCFs/chr2phase.txt
	chr_nums=($(<$chr2keep_f))

	#########################################################

	ks=($(seq 1 1 5))

	for k in ${ks[@]};
	do
		echo k=$k

		tmp1=$hic_phased_dir/generated_random_mix/all_chrs.opt.k=$k.tmp1.txt
		tmp2=$hic_phased_dir/generated_random_mix/all_chrs.opt.k=$k.tmp2.txt

		rm -f *tmp* ## so do not use parallel

		for chr_num in ${chr_nums[@]};
		do
			chr=chr$chr_num

			echo $chr

			data_chr=hic_$chr.opt.mat.txt
			shuf $data_chr | split -l $(( $(wc -l < $data_chr) * 50 / 99 )) - tmp_mat.k=$k.$chr. ## 50/99 is about 50%, and will split into two files only

			cat tmp_mat.k=$k.$chr.aa >> $tmp1
			cat tmp_mat.k=$k.$chr.ab >> $tmp2

			######################################## append pat
			data_chr=hic_$chr.opt.pat.txt
			shuf $data_chr | split -l $(( $(wc -l < $data_chr) * 50 / 99 )) - tmp_pat.k=$k.$chr. ## 50/99 is about 50%, and will split into two files only

			cat tmp_pat.k=$k.$chr.ab >> $tmp1
			cat tmp_pat.k=$k.$chr.aa >> $tmp2
		done

		for which in 1 2; do generate_hic $which $k $genome_version & done
		wait
		rm -f *tmp*		
	done
	