	#!/bin/bash
	
	wk_dir=$1
	
	###############################################
	
	## Phasing using shapeit4, b38
	## Yuanlong LIU, 30.07.2021
	
	# vcf_des_file=NA12878.Platinum.b38.vcf.gz
	# cd /mnt/ptemp/Yuanlong/unil_nas/default_nonsensitive/D2c/Yuanlong/1.Projects/15.Phasing_SV/1.Data/GM12878/mega/VCFs/
	# wget ftp://platgene_ro:''@ussd-ftp.illumina.com/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz -O $vcf_des_file
	
	# chr_nums=(X $(seq 1 1 22))
	# for chr_num in ${chr_nums[@]};
	# do
	# 	chr_name=chr$chr_num
	# 	tabix -h $vcf_des_file $chr_name | bcftools view -O b > NA12878.Platinum.b38.$chr_name.bcf
	# 	bcftools index NA12878.Platinum.b38.$chr_name.bcf
	# done
	# wait

	###############################################

	shapeit4_fun()
	{
		local wk_dir=$1		
		local chr_num=$2
		bcf2phase=$wk_dir/mega/VCFs/mega.chr$chr_num.ALL.bi.QUAL=1.hg19Tohg38.bcf
		vcf_phased=$wk_dir/mega/HapCut/shapeit4_out/mega.chr$chr_num.ALL.bi.QUAL=1.hg19Tohg38.shapedit4.vcf

		chr_name=chr$chr_num
		echo $chr_name
		genetic_map_chr=/mnt/ptemp/Yuanlong/4.Repository/1000GP_freeze9_b38/$chr_name.b38.gmap.gz
		ref_1000GP_freeze9_b38_chr=/mnt/ptemp/Yuanlong/4.Repository/1000GP_freeze9_b38/freeze9.merged.$chr_name.filtered.anno.gtonly.passonly.phased.bcf
		shapeit4 --input $bcf2phase --map $genetic_map_chr --region $chr_name --reference $ref_1000GP_freeze9_b38_chr --output $vcf_phased
	}

	###############################################

	shapeit_out=$wk_dir/mega/HapCut/shapeit4_out/
	mkdir -p $shapeit_out
	# cd $shapeit_out

	chr_nums=(X $(seq 1 1 22))
	# chr_nums=(21 22)
	for chr_num in ${chr_nums[@]};
	do
		shapeit4_fun $wk_dir $chr_num &
	done
	
  	####################################### shapeit4, https://odelaneau.github.io/shapeit4/



