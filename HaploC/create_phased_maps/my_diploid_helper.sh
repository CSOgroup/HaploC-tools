	#!/bin/bash

	## construct diploid hic map, for a single intra-chr
	## yliu
	## 16-04-2021
	wk_dir=$1
	chr_num=$2
	QUAL=$3
	max_IS=$4
	repo_dir=$5
	
 	##########################################################

    juicerDir=${repo_dir}/HaploC/dependancies/juicer/
 	
    #########################################################
	
	for chr_num_x in $chr_num; 	
	do
		echo $chr_num_x;
		chr_num_y=$chr_num_x

		merged_nodups_file=$wk_dir/mega/aligned/merged_nodup_chrXchrY/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}.snps

		chr_pos_file=$wk_dir/mega/HapCut/phased_hap_intg/hic_chr${chr_num}_QUAL=${QUAL}_maxIS=${max_IS}.chr_pos
		snp_pat_mat_file=$wk_dir/mega/HapCut/phased_hap_intg/hic_chr${chr_num}_QUAL=${QUAL}_maxIS=${max_IS}.pat_mat
		dip_out_file=$wk_dir/mega/HapCut/phased_hap_intg/hic_chr${chr_num}_QUAL=${QUAL}_maxIS=${max_IS}.dip_out

		## 10: the quality score
		prefix=$wk_dir/mega/HapCut/hic_phased/hic_chr${chr_num}_QUAL=${QUAL}_maxIS=${max_IS}

		rm -f ${dip_out_file}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.mat.txt
		rm -f ${prefix}.pat.txt
		rm -f ${prefix}.diff.txt
		rm -f ${prefix}.stat.txt

		#########################################################

		awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' ${merged_nodups_file} | $juicerDir/scripts/common/diploid.pl -s ${chr_pos_file} -o ${snp_pat_mat_file} > ${dip_out_file}
		
		mkdir -p $wk_dir/mega/HapCut/hic_phased/
		awk -v stem=$prefix. -f ${juicerDir}/scripts/common/diploid_split.awk $dip_out_file
		rm -f ${dip_out_file}

		sort -k2,2d -m ${prefix}.maternal.txt ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.mat.txt 
		sort -k2,2d -m ${prefix}.paternal.txt ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.pat.txt 

		stat_file=${prefix}.stat.txt
		mat=($(wc -l ${prefix}.maternal_both.txt))
		pat=($(wc -l ${prefix}.paternal_both.txt))
		diff=($(wc -l ${prefix}.diff.txt))

		echo "mat_both: "${mat} > ${stat_file}
		echo "pat_both: "${pat} >> ${stat_file}
		echo "diff: "${diff} >> ${stat_file}
	done
