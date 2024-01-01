	#!/bin/bash

	## keep merged_nodup.txt overlapping with SNPs
	## yliu
	## 16-04-2021
	wk_dir=$1
	chr_num_x=$2
	chr_num_y=$3	
	chr_pos_file=$4
	phased_snp_f=$5
	repo_dir=$6
	
 	##########################################################

    # juicerDir=/cso/gred/Yuanlong/3.Software_packages/juicer/
    # juicerDir=//mnt/ndata/Yuanlong/3.Packages/juicer/
 	
    #########################################################
	# for chr_num_x in $chr_num; 	
	# do
		# echo $chr_num_x;
		# chr_num_y=$chr_num_x

		merged_nodup_f=$wk_dir/aligned/merged_nodup_chrXchrY/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}
		merged_nodup_overlap_f=$wk_dir/aligned/merged_nodup_chrXchrY/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}.snps

		echo ${merged_nodup_f}
		# exit
		# merged_nodup_f=$wk_dir/aligned/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}
		# merged_nodup_overlap_f=$wk_dir/aligned/merged_nodups.txt_chr${chr_num_x}chr${chr_num_y}.snps

		# awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' ${merged_nodup_f} | $juicerDir/scripts/common/diploid.pl -s ${chr_pos_file} -o ${phased_snp_f} > ${dip_out_file}
		# awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' ${merged_nodup_f} | $juicerDir/scripts/common/diploid_keep_full_info.pl -s ${chr_pos_file} -o ${phased_snp_f} > ${merged_nodup_overlap_f}
		awk -F'[ ]' '$3>=100 && $7>=100 && $9 >= 1 && $12 >= 1' ${merged_nodup_f} | $repo_dir/HaploC/dependancies/juicer/scripts/common/diploid_keep_full_info.pl -s ${chr_pos_file} -o ${phased_snp_f} > ${merged_nodup_overlap_f}

	# done
	