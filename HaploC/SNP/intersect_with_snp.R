	#!/usr/bin/env Rscript
	args = commandArgs(trailingOnly=TRUE)
	wk_dir = args[1]
	
	chunk_num = as.numeric(args[2])
	chunk_name = sprintf('chunk%s', chunk_num)

  #################################################

	get_options = function(wk_dir)
	{
	    config_f = file.path(wk_dir, 'config/config.txt')
	    config_tab = readLines(config_f)
	    opts = as.list(sapply(strsplit(config_tab, '='), function(v) v[2]))
	    names(opts) = sapply(strsplit(config_tab, '='), function(v) v[1])
	    return(opts)
	}

	attach(get_options(wk_dir))
	source(sprintf('%s/header/HaploC_header.R', repo_dir))

	#################################################

	chr_nums2use = get_chr2use(wk_dir)
	n_chr = length(chr_nums2use)

	intra_chr_pairs = data.table(expand.grid(chr_nums2use, chr_nums2use))[Var1 == Var2]
	inter_chr_pairs = data.table(expand.grid(1:23, 1:23))[Var1 < Var2]
	inter_chr_pairs[inter_chr_pairs==23] = 'X'
	inter_chr_pairs = inter_chr_pairs[Var1 %in% chr_nums2use][Var2 %in% chr_nums2use]
  chr_pairs = rbind(intra_chr_pairs, inter_chr_pairs)
  # flag_f = file.path(wk_dir, 'mega/aligned/chunk_added.txt')

	#################################################

  run_chunk = function(chunk_name)
  {
    wk_dir_chunk = file.path(wk_dir, chunk_name)   
	  split_dir_chunk = file.path(wk_dir_chunk, 'aligned', 'merged_nodup_chrXchrY')

	  #################################################
	  
	  # unlink(split_dir_chunk, recursive = TRUE)
	  dir.create(split_dir_chunk, recursive=TRUE, showWarnings=FALSE)
		merged_nodup_file = sprintf('%s/aligned/merged_nodups.txt', wk_dir_chunk)

		################################################# split: https://www.theunixschool.com/2012/06/awk-10-examples-to-split-file-into.html
		
		## https://www.unix.com/shell-programming-and-scripting/220035-awk-script-split-file-into-multiple-files-based-many-columns.html

		setwd(split_dir_chunk)
		command_split = sprintf("awk -F'[ ]' '{print > \"merged_nodups.txt_\"$2$6 }' %s", merged_nodup_file)
		system(command_split)

		fs = list.files(pattern='merged_nodups.txt_chr*')
		for(f in fs)
		{
			split_name = strsplit(f, 'chr')[[1]]
			if(!all(split_name[-1] %in% chr_nums2use)) {unlink(f); next}						
			if(!all(split_name[-1] %in% 1:22)) next ## 'X' is in its right order
			f_new = paste0(c(split_name[1], sort(as.numeric(split_name[-1]))), collapse='chr')
			file.rename(f, f_new)
		}

		################################################# subsetting those intersecting SNPs

		no_out = foreach(i=1:nrow(intra_chr_pairs)) %do% 
		{
      chr_num = intra_chr_pairs[[1]][i]

      {
	      files = generate_diploid_mat_pat_info_n_phased(vcf_dir, chr_num)
				chr_pos_file = files['chr_pos']
				snp_pat_mat_file = files['snp_pat_mat']
      }

			command_overlap = sprintf('%s %s %s %s %s %s %s', intersect_merged_nodup, wk_dir_chunk, chr_num, chr_num, chr_pos_file, snp_pat_mat_file, repo_dir)
			cat(command_overlap, '\n')
			system(command_overlap)
			chrXchrY_file = sprintf('%s/merged_nodups.txt_chr%schr%s', split_dir_chunk, chr_num, chr_num)
			unlink(chrXchrY_file)
		}

		################################################# subsetting those intersecting SNPs

		no_out = foreach(i=1:nrow(inter_chr_pairs)) %do% 
		{
        	chr_num_x = inter_chr_pairs[[1]][i]
        	chr_num_y = inter_chr_pairs[[2]][i]

        	chr_pos_file_xy = sprintf('%s/mega.chr%s.HET.bi.QUAL=10.chr_pos', vcf_dir, c(chr_num_x, chr_num_y))
	        snp_pat_mat_file_xy = sprintf('%s/mega.chr%s.HET.bi.QUAL=10.pat_mat', vcf_dir, c(chr_num_x, chr_num_y))

        	chr_pos_file = sprintf('%s/mega.chr%schr%s.HET.bi.QUAL=10.chr_pos', vcf_dir, chr_num_x, chr_num_y)
	        snp_pat_mat_file = sprintf('%s/mega.chr%schr%s.HET.bi.QUAL=10.pat_mat', vcf_dir, chr_num_x, chr_num_y)
     

	        system(sprintf('cat %s > %s', paste0(chr_pos_file_xy, collapse=' '), chr_pos_file))
	        system(sprintf('cat %s > %s', paste0(snp_pat_mat_file_xy, collapse=' '), snp_pat_mat_file))

			command_overlap = sprintf('%s %s %s %s %s %s %s', intersect_merged_nodup, wk_dir_chunk, chr_num_x, chr_num_y, chr_pos_file, snp_pat_mat_file, repo_dir)
			cat(command_overlap, '\n')
			system(command_overlap)
			chrXchrY_file = sprintf('%s/merged_nodups.txt_chr%schr%s', split_dir_chunk, chr_num_x, chr_num_y)
			unlink(chrXchrY_file)
		}
  }
 	
 	#################################################

  vcf_dir = file.path(wk_dir, 'mega', 'VCFs')

	#################################################

  run_chunk(chunk_name)

  #################################################

  # while(1) ## previous chunks has been added; THIS strategy can be modified
  # {
  # 	Sys.sleep(sample(1:20, 1))
  # 	if(chunk_num==1 | !file.exists(flag_f)) system(sprintf('echo 0 > %s', flag_f))
  # 	chunk_added = fread(flag_f)[[1]]

  # 	if(chunk_added!=(chunk_num - 1)) next
  #   no_out = foreach(i=1:nrow(chr_pairs)) %do%
  #   {
  #   	chr_num_x = chr_pairs[[1]][i]
  #     chr_num_y = chr_pairs[[2]][i]

	#     split_dir_chunks = file.path(wk_dir, chunk_name, 'aligned', 'merged_nodup_chrXchrY')
  #   	merged_nodup_chunks = paste0(sprintf('%s/merged_nodups.txt_chr%schr%s.snps', split_dir_chunks, chr_num_x, chr_num_y), collapse=' ')
  #   	merged_nodup_mega = sprintf('%s/merged_nodups.txt_chr%schr%s.snps', split_dir_mega, chr_num_x, chr_num_y)

  #   	system(sprintf("cat %s >> %s", merged_nodup_chunks, merged_nodup_mega), wait=TRUE)

  #   	# n_lines_chunks = as.numeric( strsplit(trimws( tail(system( sprintf("wc -l %s", merged_nodup_chunks), intern=TRUE ), 1) ), ' ')[[1]][1] ) ## total lines
  #   	# n_lines_mega = as.numeric( strsplit(trimws( tail(system( sprintf("wc -l %s", merged_nodup_mega), intern=TRUE ), 1) ), ' ')[[1]][1] ) ## total lines
  #   	# if(n_lines_chunks != n_lines_mega) stop('n_lines_chunks != n_lines_mega')
  #   	# cat(i, ' done\n')
  #   }

  #   system(sprintf('echo %s > %s', chunk_num, flag_f))
  #   break
  # }

  #################################################

  wk_dir_chunk = file.path(wk_dir, chunk_name)   
 	merged_nodup_file = sprintf('%s/aligned/merged_nodups.txt', wk_dir_chunk)
 	# system(sprintf('gzip -f %s', merged_nodup_file))
 	