    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    wk_dir = args[1]

	#########################################################

	get_options = function(wk_dir)
	{
		config_f = file.path(wk_dir, 'config/config.txt')
		config_tab = readLines(config_f)
		opts = as.list(sapply(strsplit(config_tab, '='), function(v) v[2]))
		names(opts) = sapply(strsplit(config_tab, '='), function(v) v[1])
		return(opts)
	}

	opts = attach(get_options(wk_dir))
    source(sprintf('%s/header/shared_header.R', repo_dir))

	# source(sprintf('%s/header/header.R', repo_dir))
	source(sprintf('%s/header_v2.R', HaploCNV_dir))
	source(sprintf('%s/helper.R', HaploCNV_dir))
	source(sprintf('%s/HaploSharedHeader.R', HaploCNV_dir))		
	source(sprintf('%s/diffComp_helper.R', diffComp_dir))
	source(sprintf('%s/CALDER_fun.R', diffComp_dir))

	################################

	library(CALDER)
	options(scipen=25)
	registerDoParallel(cores=10)

	compartment_dir = file.path(wk_dir, "compartment")
	dir.create(file.path(compartment_dir, 'diffComp'), recursive=TRUE, showWarnings=FALSE)

	################################
	
	log_0 = try(foreach(bin_size=c(10E3, 50E3, 100E3)) %dopar% get_sparse_and_low_coverage_region(wk_dir, bin_size))

	black_list_region = try(get_black_region2(wk_dir))
	chrs2keep = unique(fread(sprintf('%s/region4comp.bed', file.path(wk_dir, 'compartment/region4comp')))$V1)
	chrs2keep = gsub('chr', '', chrs2keep)

	log_1 = try(CALDER_fun(wk_dir, map_type='phased', chr_nums2look=chrs2keep, n_cores=length(chrs2keep), black_list_region=black_list_region))
	log_2 = try(CALDER_random_fun(wk_dir, n_cores=length(chrs2keep), chr_nums2look=chrs2keep, black_list_region=black_list_region))
	log_3 = try(CALDER_fun(wk_dir, map_type='bulk', chr_nums2look=chrs2keep, n_cores=length(chrs2keep), black_list_region=NULL))
	log_4 = try(build_comp_table_opt_paired_ALL(wk_dir, chrs2keep=chrs2keep))
	log_5 = try(update_comp_Rdata(wk_dir))
	log_6 = try(diffComp_fun(wk_dir, chrs2keep))

	################################