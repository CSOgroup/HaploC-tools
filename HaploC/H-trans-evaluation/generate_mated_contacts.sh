	#!/bin/bash
	
	wk_dir=$1
	repo_dir=$2
	identity=$3

	##########################################################

    juicerDir=${repo_dir}/HaploC/dependancies/juicer/
	wk_dir_mega=$wk_dir/mega/
	merged_nodups_snp_dir=$wk_dir_mega/aligned/merged_nodup_chrXchrY

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

		prefix=$wk_dir_mega/HapCut/hic_phased/hic_chr${chr_num}_$identity

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

		######################################################### get htrans profile, both locally and globally

	    slim_out=${dip_out_f}.slim
	   
	    cut -f3,7,12,13 $dip_out_f -d" " | grep -v "mismatch" | awk 'length($3)>4 && length($4)>4' | cut -f3,4 -d" " > $slim_out

	   	sort -k2,2d -m ${prefix}.maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.mat.txt 
		sort -k2,2d -m ${prefix}.paternal_both.txt  | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > ${prefix}.pat.txt

		rm -f ${dip_out_f}
		rm -f ${prefix}.maternal.txt
		rm -f ${prefix}.paternal.txt
		rm -f ${prefix}.maternal_both.txt
		rm -f ${prefix}.paternal_both.txt
		rm -f ${prefix}.diff.txt
	}
	
	#########################################################

	for chr_num_x in ${chr_nums[@]} ; do helper_for_chr $chr_num_x & done
	wait
