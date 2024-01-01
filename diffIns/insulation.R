	#!/usr/bin/env Rscript
	args = commandArgs(trailingOnly=TRUE)
	wk_dir = args[1]
	bin_size = as.numeric(args[2])

	map_type = 'phased'

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

	################################

	registerDoParallel(cores=5)
	save_dir = file.path(wk_dir, 'diffIns')
	system(sprintf('mkdir -p %s', save_dir))

	################################
	
	insulation_score_fun = function(A, size=3, log)
	{
		insulation_score_fun_helper_log_5 = function(bin_start)
		{
			bin_mid = bin_start + size -1
			bin_end = bin_start + 2*size - 1
			up = sum(A[bin_start:bin_mid, bin_start:bin_mid])
			down = sum(A[(bin_mid+1):bin_end, (bin_mid+1):bin_end])
			all = sum(A[bin_start:bin_end, bin_start:bin_end])
			inter = 2*sum(A[bin_start:bin_mid, (bin_mid+1):bin_end])
			insulation = 1 - inter / (up + down)
			# insulation = 1 - inter / all
			# insulation = inter
			return(insulation)
		}
		# A = log2(A + 1)
		insulation_score_fun_helper = function(bin_start)
		{
			bin_mid = bin_start + size -1
			bin_end = bin_start + 2*size - 1
			up = sum(A[bin_start:bin_mid, bin_start:bin_mid])
			down = sum(A[(bin_mid+1):bin_end, (bin_mid+1):bin_end])
			all = sum(A[bin_start:bin_end, bin_start:bin_end])
			inter = 2*sum(A[bin_start:bin_mid, (bin_mid+1):bin_end])
			# insulation = 1 - inter / (up + down)
			insulation = 1 - inter / all
			# insulation = inter
			return(insulation)
		}

		insulation_score_fun_helper_log_3 = function(bin_start)
		{
			bin_mid = bin_start + size -1
			bin_end = bin_start + 2*size - 1
			up = sum(A[bin_start:bin_mid, bin_start:bin_mid])
			down = sum(A[(bin_mid+1):bin_end, (bin_mid+1):bin_end])
			all = sum(A[bin_start:bin_end, bin_start:bin_end])
			inter = 2*sum(A[bin_start:bin_mid, (bin_mid+1):bin_end])
			# insulation = 1 - inter / (up + down)
			# insulation = 1 - inter / all
			insulation = inter
			return(insulation)
		}		

		insulation_score_fun_helper_log = function(bin_start)
		{
			bin_mid = bin_start + size -1
			bin_end = bin_start + 2*size - 1
			up = sum(log(A[bin_start:bin_mid, bin_start:bin_mid]+1))
			down = sum(log(A[(bin_mid+1):bin_end, (bin_mid+1):bin_end]+1))
			inter = 2*sum(log(A[bin_start:bin_mid, (bin_mid+1):bin_end]+1))
			all = sum(log(A[bin_start:bin_end, bin_start:bin_end]+1))			
			insulation = 1 - inter / all
			# insulation = inter
			return(insulation)
		}	

		insulation_score_fun_helper_log_4 = function(bin_start)
		{
			bin_mid = bin_start + size -1
			bin_end = bin_start + 2*size - 1
			up = sum(log(A[bin_start:bin_mid, bin_start:bin_mid]+1))
			down = sum(log(A[(bin_mid+1):bin_end, (bin_mid+1):bin_end]+1))
			inter = 2*sum(log(A[bin_start:bin_mid, (bin_mid+1):bin_end]+1))
			all = sum(log(A[bin_start:bin_end, bin_start:bin_end]+1))			
			# insulation = 1 - inter / all
			insulation = inter
			return(insulation)
		}				

		bin_starts = 1:(nrow(A) - 2*size + 1)
		if(log==0) insulations = sapply(bin_starts, function(bin_start) insulation_score_fun_helper(bin_start))
		if(log==1) insulations = sapply(bin_starts, function(bin_start) insulation_score_fun_helper_log(bin_start))		
		if(log==3) insulations = sapply(bin_starts, function(bin_start) insulation_score_fun_helper_log_3(bin_start))		
		if(log==4) insulations = sapply(bin_starts, function(bin_start) insulation_score_fun_helper_log_4(bin_start))		
		if(log==5) insulations = sapply(bin_starts, function(bin_start) insulation_score_fun_helper_log_5(bin_start))		

		# insulations = foreach(bin=bin_starts, .combine='c') %dopar% insulation_score_fun_helper(bin)
		# insulations = as.vector(scale(insulations))
		names(insulations) = rownames(A)[size:(nrow(A) - size)]
		insulation_tab = data.table(bin_index=as.numeric(names(insulations)), ins=insulations)
		return(insulation_tab)
	}

	ins_core = function(hic_f, ins_bed)
	{
		insulation_bed = foreach(chr_num_x=chr_nums2use, .combine=rbind) %dopar%
		{
			chr_name = paste0("chr", chr_num_x)
			chr2straw = as.character(chr_num_x)
			if(any(grepl('chr', readHicChroms(hic_f)$name))) chr2straw = chr_name

			dump_mat = try(silent=TRUE, straw(matrix='oe', binsize=bin_size, chr1loc=chr2straw, chr2loc=chr2straw, unit='BP', fname=hic_f, norm='KR'))
			if(class(dump_mat)=='try-error') dump_mat = try(silent=TRUE, straw(matrix='oe', binsize=bin_size, chr1loc=chr2straw, chr2loc=chr2straw, unit='BP', fname=hic_f, norm='VC_SQRT'))
			
			if(class(dump_mat)[1]=='try-error')
			{
				insulation_bed = data.table(chr=chr_name, start=0, end=0, score=0)
				return(insulation_bed)
			}

			if(nrow(na.omit(dump_mat)) < 100)
			{
				insulation_bed = data.table(chr=chr_name, start=0, end=0, score=0)
				return(insulation_bed)
			}

			dump_mat = as.data.table(dump_mat)

	        colnames(dump_mat) = c("pos_1", "pos_2", "counts")
			dump_mat = dump_mat[!is.nan(counts)]

	        mat = convert_contact_format_fun(dump_mat, from='format_dump', to='format_sparse', chr_num=chr_num_x, binsize=bin_size)$format_sparse
	        mat[abs(row(mat) - col(mat)) > 20] = 0 ## 50 bins
	        n_row = nrow(mat)
	        insulation_tab = insulation_score_fun(mat, log=0, size=20) 
	        # print('WARNING: seems size=10 was used before June 2023! CHECK')
	        start = as.character(bin_size*(insulation_tab$bin_index - 1) + 1)
	        end = as.character(bin_size*(insulation_tab$bin_index))
	        # insulation_bed = data.table(chr=chr_name, start=start, end=end, score=insulation_tab$ins - 0.5)
	        insulation_bed = data.table(chr=chr_name, start=start, end=end, score=insulation_tab$ins - 0.5)	        
	        insulation_bed
		}
		
		insulation_bed = insulation_bed[!is.nan(score)]
        fwrite(insulation_bed, file=ins_bed, sep='\t', quote=FALSE, col.names=FALSE)
        return(insulation_bed)
	}

	################################

	if(is.null(bin_size) | is.na(bin_size)) bin_size = 40E3

	################################

	if(map_type=="bulk") identities = "bulk"
	if(map_type=="phased" | map_type=='rand') identities = mates

	chr_nums2use = get_chr2use(wk_dir)

	################################

	ins_mates = foreach(mate=identities) %do%
	{
		if(map_type!='rand')
		{
			ins_bed = sprintf("%s/ins.%s.bedgraph", save_dir, mate)
			hic_f = get_mate_hic(wk_dir, mate)
			if(map_type=='bulk') hic_f = get_bulk_hic(wk_dir)
			ins_tab = ins_core(hic_f, ins_bed)
			colnames(ins_tab)[4] = mate
			return(ins_tab)
		}


		if(map_type=='rand') 
		{
			foreach(k=1:5) %do%
			{
				cat('k=', k, '\n')
				ins_bed = sprintf('%s/rand_k=%s.ins.%s.bedgraph', ins_rand_dir, k, mate)
				hic_f = sprintf('%s/mega/HapCut/hic_phased/generated_random_mix/all_chrs.opt.k=%s.rand%s.hic', wk_dir, k, ifelse(mate=='mat', 1, 2))
				try(ins_core(hic_f, ins_bed))
			}
		}       
	}

	if(map_type=='phased')
	{	
		ins_mates_tab = Reduce(merge, ins_mates)
		ins_mates_tab = ins_mates_tab[mat > 0 & pat > 0][order(chr, as.numeric(start))]
		ins_mates_tab$diffIns = ins_mates_tab$pat - ins_mates_tab$mat
		ins_mates_tab$log2diffIns = log2(ins_mates_tab$pat/ins_mates_tab$mat)


		# Apply a moving average filter / chatGPT
		window_size <- 5
		smoothed_signal <- filter(ins_mates_tab$diffIns, rep(1/window_size, window_size), sides = 2)
		ins_mates_tab$diffIns = smoothed_signal

		diffIns_bed = sprintf("%s/diffIns.bedgraph", save_dir)
		log2diffIns_bed = sprintf("%s/log2diffIns.bedgraph", save_dir)

		fwrite(ins_mates_tab[, c(1:3, 6)], file=diffIns_bed, sep='\t', quote=FALSE, col.names=FALSE)
		fwrite(ins_mates_tab[, c(1:3, 7)], file=log2diffIns_bed, sep='\t', quote=FALSE, col.names=FALSE)		
	}

	################################

	log2show = sprintf("Done diffIns for %s", basename(wk_dir))
	cat(log2show, "\n")
	