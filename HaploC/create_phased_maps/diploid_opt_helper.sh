	#!/bin/bash

	## construct diploid hic map, for a single intra-chr
	## yliu
	## 16-04-2021
	wk_dir=$1
	# chr_nums=$1
	genome_version=$2
 	##########################################################
 	
 	# wk_dir=/mnt/ptemp/Yuanlong/1.Data/phasing/IMR90/HiC050_056/mega/
    juicerDir=/mnt/ndata/Yuanlong/3.Packages/juicer/
	juicer_jar=/mnt/etemp/Yuanlong/nas_yliu/1.HiC/12.Data/1.Test_juice_pre/juicer_tools_1.19.02.jar

    #########################################################
    identity="ground_truth"
    # identify="opt"

	for chr_num_x in {1..22} X;
	# for chr_num_x in $chr_nums;	
	do
		echo $chr_num_x;
		chr_num_y=$chr_num_x
		chr_num=$chr_num_x

		merged_nodups_file=$wk_dir/aligned/merged_nodup_chrXchrY/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}.snps

		chr_pos_file=$wk_dir/HapCut/phased_hap_intg/hic_chr${chr_num}_$identity.chr_pos
		snp_pat_mat_file=$wk_dir/HapCut/phased_hap_intg/hic_chr${chr_num}_$identity.pat_mat
		dip_out_file=$wk_dir/HapCut/phased_hap_intg/hic_chr${chr_num}_$identity.dip_out

		prefix=$wk_dir/HapCut/hic_phased/hic_chr${chr_num}_$identity
		mkdir -p $wk_dir/HapCut/hic_phased/

		rm -f ${dip_out_file}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.mat.txt
		rm -f ${prefix}.pat.txt
		rm -f ${prefix}.diff.txt

		echo $snp_pat_mat_file

		awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' $merged_nodups_file | $juicerDir/scripts/common/diploid.pl -s $chr_pos_file -o $snp_pat_mat_file > $dip_out_file
		
		awk -v stem=$prefix. -f ${juicerDir}/scripts/common/diploid_split.awk $dip_out_file
		rm $dip_out_file
		# awk -v stem=$prefix. -f ${juicerDir}/scripts/common/diploid_split_inter_YL.awk $dip_out_file

		sort -k2,2d -m ${prefix}.maternal.txt ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.mat.txt 
		sort -k2,2d -m ${prefix}.paternal.txt ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.pat.txt 

		# java -Xmx200g -jar $juicer_jar pre -r 20000,25000,40000,50000,100000,500000,1000000 ${prefix}.mat.txt ${prefix}.mat.hic hg19
		# java -Xmx200g -jar $juicer_jar pre -r 20000,25000,40000,50000,100000,500000,1000000 ${prefix}.pat.txt ${prefix}.pat.hic hg19
	done

	generate_hic()
	{
		local wk_dir=$1
		local which=$2
		cat $wk_dir/HapCut/hic_phased/hic_chr*$identity.$which.txt > $wk_dir/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.txt
		java -Xmx200g -jar $juicer_jar pre -r 10000,20000,25000,40000,50000,100000,500000,1000000 $wk_dir/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.txt $wk_dir/HapCut/hic_phased/hic_ALL_chrs.$identity.$which.hic $genome_version
	}

	for which in mat pat; do generate_hic $wk_dir $which & done

	


