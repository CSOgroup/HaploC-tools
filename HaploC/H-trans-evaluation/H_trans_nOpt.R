	# !/usr/bin/env Rscript

	args = commandArgs(trailingOnly=TRUE)
	wk_dir = args[1]

	repo_dir = '/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/0.Scripts/'    
	source(sprintf('%s/header/HaploC_header.R', repo_dir))

	####################################
	
    chr_nums2use = get_chr2use(wk_dir)
    registerDoParallel(cores=length(chr_nums2use))

    ####################################

	QUAL = 10; max_IS = '50'
	identity = sprintf('QUAL=10_maxIS=%s', max_IS)

	####################################

	block_ratio_thresh = 0.01 ## it was 0.05 before 02.08.2021

	no_out = foreach(chr_num=chr_nums2use) %dopar%
	{
		vcf_file = sprintf('%s/mega/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.hap.phased.VCF', wk_dir, chr_num, QUAL, max_IS)
		vcf_tab = try(vcf2tab(block_ratio_thresh=block_ratio_thresh, vcf_file, which_blocks=NA, gz=FALSE, phased=TRUE, save2wig=TRUE, consider_ps=TRUE, wig_ground_truth=NA, format_shapeit2_hap=FALSE))
		keep_pos = vcf_tab$POS
		snps_used_1 = generate_diploid_mat_pat_info_phased_vcf_tab(vcf_tab, file.path(wk_dir, 'mega'), chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA, keep_pos=keep_pos, keep_range=NA)

		command_diploid = sprintf('%s %s %s %s %s %s', my_diploid_helper, wk_dir, chr_num, QUAL, max_IS, repo_dir)
		system(command_diploid, wait=TRUE)

		cmd_htrans = sprintf('%s %s %s %s %s generate_data', Rscript, htrans_cis_info_global_fun, wk_dir, identity, chr_num)
		system(cmd_htrans)
	}


    #################################### also do for maxIS=20Mb

	QUAL = 10; max_IS = '20'
	identity = sprintf('QUAL=10_maxIS=%s', max_IS)

	####################################

	block_ratio_thresh = 0.01 ## it was 0.05 before 02.08.2021

	no_out = foreach(chr_num=chr_nums2use) %dopar%
	{
		vcf_file = sprintf('%s/mega/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.hap.phased.VCF', wk_dir, chr_num, QUAL, max_IS)
		vcf_tab = try(vcf2tab(block_ratio_thresh=block_ratio_thresh, vcf_file, which_blocks=NA, gz=FALSE, phased=TRUE, save2wig=TRUE, consider_ps=TRUE, wig_ground_truth=NA, format_shapeit2_hap=FALSE))
		keep_pos = vcf_tab$POS
		snps_used_1 = generate_diploid_mat_pat_info_phased_vcf_tab(vcf_tab, file.path(wk_dir, 'mega'), chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA, keep_pos=keep_pos, keep_range=NA)

		command_diploid = sprintf('%s %s %s %s %s %s', my_diploid_helper, wk_dir, chr_num, QUAL, max_IS, repo_dir)
		system(command_diploid, wait=TRUE)

		cmd_htrans = sprintf('%s %s %s %s %s generate_data', Rscript, htrans_cis_info_global_fun, wk_dir, identity, chr_num)
		system(cmd_htrans)
	}

	