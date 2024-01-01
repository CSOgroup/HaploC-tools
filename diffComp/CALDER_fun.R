
	################################ update opt mated opt comp

	build_comp_table_opt_paired_ALL = function(wk_dir, bin_sizes=c(50E3, 100E3, 150E3, 200E3), with_ref=TRUE, chrs2keep)
	{
		IDs = c('phased', sprintf("random_k=%s", 1:5))
		for(ID in IDs) build_comp_table_opt_paired(wk_dir, ID, chrs2keep=chrs2keep)
	}

	build_comp_table_opt_paired = function(wk_dir, ID, bin_sizes=c(50E3, 100E3, 150E3, 200E3), with_ref=TRUE, chrs2keep) 
	{
		get_opt_index_mates = function(consist_tab_mates)
		{
			n = nrow(consist_tab_mates[['mat']])
			threshs = c(0.05, 0.1, 0.15)
			candidate_indices_tab = foreach(thresh=threshs, .combine=cbind) %do%
			{
				candidate_indices = foreach(index=1:n, .combine=c) %do%
				{
					tmp = rbind(consist_tab_mates$mat[index, ], consist_tab_mates$pat[index, ])
					candidate_tab = 1*(tmp[, -1] > apply(tmp[, -1], 1, max) - thresh)
					candidate = min(which(candidate_tab[1, ] * candidate_tab[2, ] == 1))
					candidate
				}
				candidate_indices
			}

			colnames(candidate_indices_tab) = threshs
			opt_index = apply(candidate_indices_tab, 1, min)
			return(opt_index)
		}

		################################

		consist_tab_mates = foreach::foreach(mate=mates) %do%
		{
			CALDER_dir = file.path(wk_dir, 'compartment', ID, mate)
			consist_tab = foreach::foreach(bin_size=bin_sizes, .combine=merge) %do%
			{
				bin_size_kb = sprintf("%skb", bin_size/1E3)
        		CALDER_dir_binsize = file.path(CALDER_dir, 'intermediate_data/sub_compartments', bin_size_kb)
				consist_tab_tmp = data.table::data.table(chr=chrs2keep, val=0)
				consist_tab_tmp$val = foreach::foreach(chr=chrs2keep, .combine=c) %do%
				{
				    log_file = sprintf('%s/chr%s_log.txt', CALDER_dir_binsize, chr)
				    cor_val = as.numeric(strsplit(readLines(log_file)[5], 'this chr is:')[[1]][2])
				    cor_val
				}
				colnames(consist_tab_tmp)[2] = bin_size_kb
				consist_tab_tmp
			}

			consist_tab = consist_tab[order(match(consist_tab$chr, chrs2keep)), ]
			consist_tab
		}

		names(consist_tab_mates) = mates
		opt_index = get_opt_index_mates(consist_tab_mates)

		################################ Choose the best bin size (as small as possible, and save the chosen comp to the opt dir)

		silent_out = foreach::foreach(mate=mates) %do% ## do not use dopar
		{
			CALDER_dir = file.path(wk_dir, 'compartment', ID, mate)
			CALDER_dir_opt = file.path(CALDER_dir, "sub_compartments")
			CALDER_dir_opt_n_mate = paste0(CALDER_dir, '.n_mate')

			if(file.exists(CALDER_dir_opt_n_mate)) system(sprintf('rm -r %s', CALDER_dir_opt_n_mate))
			system(sprintf('cp -r %s %s', CALDER_dir, CALDER_dir_opt_n_mate))
			system(sprintf('rm -r %s/intermediate_data', CALDER_dir_opt_n_mate))
			
			consist_tab = consist_tab_mates[[mate]]		

			names(opt_index) = gsub(":", "", consist_tab$chr)
			opt_bin_tab = data.table::data.table(chr=names(opt_index),  opt_binsize=(colnames(consist_tab)[-1])[opt_index])

			consist_tab_tmp = cbind( consist_tab, opt_bin_tab[, 'opt_binsize'])
			colnames(consist_tab_tmp)[ncol(consist_tab_tmp)] = 'opt_binsize'

			cor_all_file = file.path(CALDER_dir_opt, "cor_with_ref.ALL.txt")
			data.table::fwrite(consist_tab_tmp, file=cor_all_file, col.names=TRUE, sep="\t")

			cor_with_ref = foreach::foreach(j=1:nrow(opt_bin_tab), .combine=c) %do%
			{
				chr_query = opt_bin_tab$chr[j]
				opt_binsize = opt_bin_tab$opt_binsize[j]
				if(!grepl('kb', opt_binsize)) return(0)

				row_index = which(consist_tab$chr == chr_query)
				col_index = which(colnames(consist_tab)==opt_binsize)
				as.data.frame(consist_tab)[row_index, col_index] ## some wired behavior when using data.table. Convert to data.frame
			}

			cor_with_ref_tab = data.table::data.table(chr=consist_tab$chr, cor=cor_with_ref, chr_num=gsub("chr", "", opt_bin_tab$chr))
			cor_with_ref_tab = cor_with_ref_tab[order(match(consist_tab$chr, chrs2keep)), ]


			cor_opt_file = file.path(CALDER_dir_opt, "cor_with_ref.txt")
			data.table::fwrite(cor_with_ref_tab[, 1:2], file=cor_opt_file, col.names=TRUE, sep="\t")

			################################ make plot

			cor_with_ref_tab$index = 1:nrow(cor_with_ref_tab)
			cor_opt_plot_file = file.path(CALDER_dir_opt, "cor_with_ref.pdf")


		    pdf(cor_opt_plot_file, width=8, height=6)
		    p = ggplot2::ggplot(data=cor_with_ref_tab, ggplot2::aes(x=index, y=cor, label = format(round(cor_with_ref_tab$cor, 2), nsmall = 2))) + ggplot2::geom_line(color="red") + ggplot2::geom_point()
		    p = p + ggplot2::geom_label() + ggplot2::scale_x_continuous(labels=cor_with_ref_tab$chr_num, breaks=1:nrow(cor_with_ref_tab))
		    p = p + ggplot2::geom_hline(yintercept=ifelse(with_ref, 0.4, 0.15), linetype=2, color='darkblue')
		    if(with_ref) p = p + ggplot2::xlab('chr') + ggplot2::ylab('Rho') + ggplot2::ggtitle('Correlation of compartment rank with reference\nRho < 0.4 indicates imprecise compartment calling')
		    if(!with_ref) p = p + ggplot2::xlab('chr') + ggplot2::ylab('Rho') + ggplot2::ggtitle('Correlation of compartment rank with reference\nRho < 0.15 indicates imprecise compartment calling')
		    print(p)    
		    dev.off()			

		    ################################

			opt_sub_comp_beds = foreach::foreach(i=1:nrow(opt_bin_tab), .combine=c) %do%
			{
				opt_bin_size_kb = opt_bin_tab[[2]][i]
				if(is.na(opt_bin_size_kb)) return(NULL)
				chr_num = opt_bin_tab[[1]][i]
				opt_sub_comp_bed_chr = sprintf("%s/chr%s_sub_compartments.bed", file.path(CALDER_dir, 'intermediate_data/sub_compartments', opt_bin_size_kb), chr_num)
				opt_sub_comp_bed_chr
			}

			# save_opt_file = file.path(CALDER_dir_opt, "all_sub_compartments.bed")
			save_opt_file = file.path(CALDER_dir_opt, sprintf("all_sub_compartments.bed"))

			cmd_opt = sprintf("cat %s > %s", paste0(opt_sub_comp_beds, collapse=" "), save_opt_file)
			system(cmd_opt)


			comp_raw = data.table::fread(save_opt_file)
			comp_tab = make_comp_tab(comp_raw, bin_size=10E3)
            cols2keep = c("chr", "pos_start", "pos_end", "comp_name", "comp_rank", "continous_rank")
            comp_tab = comp_tab[, c("chr", "pos_start", "pos_end", "comp_name", "comp_rank", "continous_rank")]
            save_opt_bed_file = file.path(CALDER_dir_opt, sprintf("all_compartment_rank.bedgraph"))
            comp_tab$continous_rank = comp_tab$continous_rank - 0.5
            write.table(comp_tab[, c("chr", "pos_start", "pos_end", "continous_rank")], file = save_opt_bed_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

			save_opt_tab_file = file.path(CALDER_dir_opt, sprintf("all_sub_compartments.tsv"))
			data.table::fwrite(comp_tab, file=save_opt_tab_file, sep='\t')
		}

		rank_tab = foreach(mate=mates, .combine=merge) %do%
		{
			CALDER_dir = file.path(wk_dir, 'compartment', ID, mate)
			old_bed = file.path(CALDER_dir, 'sub_compartments/all_sub_compartments.bed')
			new_bed = file.path(CALDER_dir, 'all_sub_compartments.%s.bed')
			old_bedgraph = file.path(CALDER_dir, 'sub_compartments/all_compartment_rank.bedgraph')
			new_bedgraph = file.path(CALDER_dir, sprintf('all_compartment_rank.%s.bedgraph', mate))	
			
			system(sprintf('mv %s %s', old_bed, new_bed))	
			system(sprintf('mv %s %s', old_bedgraph, new_bedgraph))	

			cmd_replaceX = sprintf("sed -i 's/chrNA/chrX/g' %s", new_bed) ## the bed file contains NA because of chrX
			system(try(cmd_replaceX))
			cmd_replaceX = sprintf("sed -i 's/chrNA/chrX/g' %s", new_bedgraph) ## the bed file contains NA because of chrX
			system(try(cmd_replaceX))

			rank_tab = fread(new_bedgraph)	
			rank_tab = expand_comp_tab(rank_tab)[, c(1,6,7,4)]
			colnames(rank_tab)[4] = mate
			rank_tab
		}

		{
			rank_tab$diff = rank_tab$pat - rank_tab$mat
			diff_bedgraph = file.path(wk_dir, 'compartment/diffComp', sprintf('%s.diff_rank.bedgraph', ID))
			unlink(diff_bedgraph)
			write.table(rank_tab[, c(1,2,3,6)], file=diff_bedgraph, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

			# diff_bedgraph_n_rand = file.path(wk_dir, 'CALDER', sprintf('%s.diff_rank.bedgraph', ID))
			# write.table(rank_tab[, c(1,2,3,6)], file=diff_bedgraph_n_rand, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')			
		}

	}

	################################

	CALDER_fun = function(wk_dir, map_type, chr_nums2look, n_cores, black_list_region)
	{
		if(map_type=="bulk") identities = "bulk"
		if(map_type=="phased") identities = mates

		bin_sizes = c(50E3, 100E3, 150E3, 200E3)
		bin_size = bin_sizes[1]

		################################

		rank_tab = foreach(mate=identities, .combine=merge) %do%
		{
			hic_f = get_mate_hic(wk_dir, mate)			
			if(mate=='bulk') hic_f = get_bulk_hic(wk_dir)

			if(!file.exists(hic_f)) stop('not exist')	
			CALDER_dir = file.path(wk_dir, 'compartment', map_type, mate)
			if(mate=='bulk') CALDER_dir = file.path(wk_dir, 'compartment', mate)			
			unlink(CALDER_dir, recursive=TRUE)

			CALDER(contact_file_hic=hic_f, chrs=chr_nums2look, bin_size=bin_size, bin_sizes=bin_sizes,  genome='hg19', save_dir=CALDER_dir, save_intermediate_data=FALSE, n_cores=n_cores, sub_domains=FALSE, black_list_region=black_list_region)
			
			old_bed = file.path(CALDER_dir, 'sub_compartments/all_sub_compartments.bed')
			new_bed = file.path(CALDER_dir, sprintf('all_sub_compartments.%s.bed', mate))

			old_bedgraph = file.path(CALDER_dir, 'sub_compartments/all_compartment_rank.bedgraph')
			new_bedgraph = file.path(CALDER_dir, sprintf('all_compartment_rank.%s.bedgraph', mate))	
			
			system(sprintf('mv %s %s', old_bed, new_bed))	
			system(sprintf('mv %s %s', old_bedgraph, new_bedgraph))	

			cmd_replaceX = sprintf("sed -i 's/chrNA/chrX/g' %s", new_bed) ## the bed file contains NA because of chrX
			system(try(cmd_replaceX))
			cmd_replaceX = sprintf("sed -i 's/chrNA/chrX/g' %s", new_bedgraph) ## the bed file contains NA because of chrX
			system(try(cmd_replaceX))

			rank_tab = fread(new_bedgraph)	
			rank_tab = expand_comp_tab(rank_tab)[, c(1,6,7,4)]
			colnames(rank_tab)[4] = mate
			rank_tab
		}

		if(map_type=='phased')
		{
			rank_tab$diff = rank_tab$pat - rank_tab$mat
			# diff_bedgraph = file.path(wk_dir, 'CALDER/diffComp/diff_rank.bedgraph')
			unlink(diff_bedgraph)
			write.table(rank_tab[, c(1,2,3,6)], file=diff_bedgraph, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		}
	}

	###########################################

	CALDER_random_fun = function(wk_dir, chr_nums2look, n_cores, black_list_region)
	{
		no_out = foreach(k=1:5) %do%
		{
			line = sprintf("random_k=%s", k)

			################################

			rank_tab = foreach(mate=mates, .combine=merge) %do%
			{
				if(mate=='mat') hic_f = sprintf('%s/mega/HapCut/hic_phased/generated_random_mix/all_chrs.opt.k=%s.rand1.hic', wk_dir, k)
				if(mate=='pat') hic_f = sprintf('%s/mega/HapCut/hic_phased/generated_random_mix/all_chrs.opt.k=%s.rand2.hic', wk_dir, k)
				if(!file.exists(hic_f)) return('not exist')		

				CALDER_dir = file.path(wk_dir, 'compartment', line, mate)

				CALDER(contact_file_hic=hic_f, chrs=chr_nums2look, bin_size=50E3, bin_sizes=c(50E3, 100E3, 150E3, 200E3, 250E3),  genome='hg19', save_dir=CALDER_dir, save_intermediate_data=FALSE, n_cores=n_cores, sub_domains=FALSE, black_list_region=black_list_region)
				old_bed = file.path(CALDER_dir, 'sub_compartments/all_sub_compartments.bed')
				new_bed = file.path(CALDER_dir, sprintf('all_sub_compartments.%s.bed', mate))

				old_bedgraph = file.path(CALDER_dir, 'sub_compartments/all_compartment_rank.bedgraph')
				new_bedgraph = file.path(CALDER_dir, sprintf('all_compartment_rank.%s.bedgraph', mate))	
				
				system(sprintf('mv %s %s', old_bed, new_bed))	
				system(sprintf('mv %s %s', old_bedgraph, new_bedgraph))	

				cmd_replaceX = sprintf("sed -i 's/chrNA/chrX/g' %s", new_bed) ## the bed file contains NA because of chrX
				system(try(cmd_replaceX))
				cmd_replaceX = sprintf("sed -i 's/chrNA/chrX/g' %s", new_bedgraph) ## the bed file contains NA because of chrX
				system(try(cmd_replaceX))

				rank_tab = fread(new_bedgraph)	
				rank_tab = expand_comp_tab(rank_tab)[, c(1,6,7,4)]
				colnames(rank_tab)[4] = mate
				rank_tab

			}

			rank_tab$diff = rank_tab$pat - rank_tab$mat
			diff_bedgraph = file.path(wk_dir, sprintf('compartment/diffComp/%s.diff_rank.bedgraph', line))
			write.table(rank_tab[, c(1,2,3,6)], file=diff_bedgraph, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		}

	}

	################################

	update_comp_Rdata = function(wk_dir)
	{
		lines = c('phased', sprintf("random_k=%s", 1:5)) ## do not include line without random
		comp_line_f = file.path(wk_dir, "/compartment/comp_tab_line.Rdata")
		comp_line_cor_f = file.path(wk_dir, "/compartment/comp_tab_line.cor.Rdata")

		comp_tab_line = comp_tab_line_cor = list()


		map_type = 'bulk'
		CALDER_dir = file.path(wk_dir, 'compartment', map_type)
		opt_comp_f = file.path(CALDER_dir, sprintf('all_sub_compartments.%s.bed', map_type))
		
		cor_with_ref_f = file.path(CALDER_dir, '/sub_compartments/cor_with_ref.txt')
		cor_with_ref_tab = fread(cor_with_ref_f)
		cor_with_ref_tab = cor_with_ref_tab[!grepl('chrX', chr)]
		cor_with_ref_tab$chr_num = as.numeric(gsub('chr', '', cor_with_ref_tab$chr))
		cor_with_ref_tab$chr = paste0(cor_with_ref_tab$chr, ':')

		comp_raw = fread(opt_comp_f)
		comp_tab = make_comp_tab(comp_raw, bin_size=10E3)
		comp_tab_line[[map_type]] = comp_tab
		comp_tab_line_cor[[map_type]] = cor_with_ref_tab


		for(line in lines)
		{
			no_out = foreach(mate=c(mates)) %do% ## do not use dopar
			{
				CALDER_dir = file.path(wk_dir, 'compartment', line, mate)
				opt_comp_f = file.path(CALDER_dir, sprintf('all_sub_compartments.%s.bed', mate))
				
				cor_with_ref_f = file.path(CALDER_dir, sprintf('/sub_compartments/cor_with_ref.txt', mate))
				cor_with_ref_tab = fread(cor_with_ref_f)
				cor_with_ref_tab = cor_with_ref_tab[!grepl('chrX', chr)]
				cor_with_ref_tab$chr_num = as.numeric(gsub('chr', '', cor_with_ref_tab$chr))
				cor_with_ref_tab$chr = paste0(cor_with_ref_tab$chr, ':')

				comp_raw = fread(opt_comp_f)
				comp_tab = make_comp_tab(comp_raw, bin_size=10E3)
				name = sprintf("%s.%s", line, mate)
				comp_tab_line[[name]] = comp_tab
				comp_tab_line_cor[[name]] = cor_with_ref_tab
			}
		}

		save(comp_tab_line, file=comp_line_f)
		save(comp_tab_line_cor, file=comp_line_cor_f)
	}
	