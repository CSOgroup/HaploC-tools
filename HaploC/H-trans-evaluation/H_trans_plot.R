	# !/usr/bin/env Rscript

	args = commandArgs(trailingOnly=TRUE)
	wk_dir = args[1]

	repo_dir = '/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/0.Scripts/'    
	source(sprintf('%s/header/HaploC_header.R', repo_dir))

	##################################### generate H-trans for plot

	identity = c('QUAL=10_maxIS=20', 'opt')[1]
	cmd_htrans_plot = sprintf('%s %s %s %s psudo_chr generate_plot', Rscript, htrans_cis_info_global_fun, wk_dir, identity)
	system(cmd_htrans_plot)

	identity = c('QUAL=10_maxIS=50', 'opt')[1]
	cmd_htrans_plot = sprintf('%s %s %s %s psudo_chr generate_plot', Rscript, htrans_cis_info_global_fun, wk_dir, identity)
	system(cmd_htrans_plot)

	#####################################

	