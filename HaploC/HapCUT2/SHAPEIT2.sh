	#!/bin/bash
	
	chr=$1
	vcf_chr_char=$2
	frag_hapcut_chr=$3
	frag_intg_chr=$4
	shapeit_out=$5
	append=$6
	repo_dir=$7

  	####################################### shapeit2

  	shapeit=$repo_dir/HaploC/dependancies/SHAPEIT2/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit
	thousandG3_dir=$repo_dir/genomicData/1000GP_Phase3/
	mkdir -p ${shapeit_out}
	cd ${shapeit_out}

	samplehaps=$repo_dir/HaploC/dependancies/IntegratedPhasing/samplehaps.py
	encodereads=$repo_dir/HaploC/dependancies/IntegratedPhasing/encodereads.py

	####################################### SHOULD THIS PART ALSO BE X-SPECIFIC? NON-PAR? 05-Aug-2023 / changed hap and legend to NONPAR
	
	echo chr${chr}
	reference_panel_hap=$thousandG3_dir/1000GP_Phase3_chr${chr}.hap.gz
	reference_panel_legend=$thousandG3_dir/1000GP_Phase3_chr${chr}.legend.gz ## 
	reference_panel_samples=$thousandG3_dir/1000GP_Phase3.sample ## HG00096 GBR EUR male

	if [[ $chr == "X" ]]
	then
	   	genetic_map=$thousandG3_dir/genetic_map_chrX_nonPAR_combined_b37.txt ## 154969627 0.658886078468027 0.000396779120536043
		reference_panel_hap=$thousandG3_dir/1000GP_Phase3_chrX_NONPAR.hap.gz
		reference_panel_legend=$thousandG3_dir/1000GP_Phase3_chrX_NONPAR.legend.gz
	else
	   genetic_map=$thousandG3_dir/genetic_map_chr${chr}_combined_b37.txt ## 154969627 0.658886078468027 0.000396779120536043
	fi
	echo $genetic_map

	#######################################

	vcf_base_name=`basename "$vcf_chr_char"`

	if [[ "$append" == 'create' ]]; then

		echo "Will do: generating pseudo_reads"
		rand=$RANDOM
		prefix=${vcf_base_name}_delete_${rand}_`date '+%F_%H:%M:%S'`
		vcf_chr_num=${vcf_chr_char}.${prefix}

		echo ${vcf_base_name}
		sed "s/chr${chr}/${chr}/g" $vcf_chr_char > $vcf_chr_num

		$shapeit -check --input-vcf $vcf_chr_num -R $reference_panel_hap $reference_panel_legend $reference_panel_samples --output-log $prefix
		$shapeit --input-vcf $vcf_chr_num -R $reference_panel_hap $reference_panel_legend $reference_panel_samples -M $genetic_map --output-graph ${prefix}.graph --seed 123456789 --exclude-snp ${prefix}.snp.strand.exclude
		$shapeit -convert --input-graph ${prefix}.graph --output-max ${prefix}
		mv ${prefix}.haps ${vcf_base_name}.haps

		python2.7 $samplehaps $vcf_chr_num $prefix 1000 $repo_dir
		python2.7 $encodereads $prefix.hapsamples > ${vcf_base_name}.pseudo_reads
		rm ${vcf_chr_num}
		rm *${prefix}*
		echo "Finished: creating pseudo_reads"			
	else
		echo "Will do: appending pseudo_reads"
		cp ${frag_hapcut_chr} ${frag_intg_chr}

		vcf_base_name_ALL=${vcf_base_name/HET/ALL} ## because the pseudo_reads is named as "ALL"
		cat ${vcf_base_name_ALL}.pseudo_reads >> ${frag_intg_chr}

		# cat ${vcf_base_name}.pseudo_reads >> ${frag_intg_chr}
		echo "Finished: appending pseudo_reads"
	fi
