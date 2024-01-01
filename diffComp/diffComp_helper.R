    get_black_region2 = function(wk_dir) ## deal with regions with translocation / low inter-region contacts
    {
        ref = 'hg19'
        bin_size = 50E3 ## not sure about this
        CNV_dir = file.path(wk_dir, sprintf("CNV/%skb/", bin_size/1E3))

        chr_nums2use = get_chr2use(wk_dir)

        green_list_f =  green_region_bed = file.path(CNV_dir, sprintf("green_list.bed"))
        green_list = fread(green_list_f)
        colnames(green_list) = GR_3_cols

        region4comp = get_largest_block(wk_dir=wk_dir, green_list=green_list, chr_nums=chr_nums2use, thresh_mate=0.05) ## by keeping largest block
        region4comp = fread(sprintf('%s/region4comp.bed', file.path(wk_dir, 'compartment/region4comp')))
        colnames(region4comp) = GR_3_cols
        region4compGR = makeGRangesFromDataFrame(region4comp)

        chr_size = fread(chrom_sizes)
        colnames(chr_size) = c("chr", "len")
        chr_size = chr_size[chr %in% paste0("chr", c(1:22, "X"))]
        chr_size_GR = makeGRangesFromDataFrame(data.table(chr=chr_size$chr, start=1, end=chr_size$len))
        black_region = setdiff(chr_size_GR, region4compGR)
        black_region = as.data.table(black_region)[, 1:3]
        colnames(black_region) = GR_3_cols
        black_region$end = black_region$end + 1
        black_region$chr = as.character(black_region$chr)

        return(black_region)
    }

    get_largest_block = function(bin_size=100E3, wk_dir, green_list=NULL, chr_nums, thresh_mate) ## if density < 0.2 / largest block for deriving compartment / no much need to correct for CNV, because the criteria seems powerful to handle CNV / 25-02-2023
    {

        region4comp_dir = file.path(wk_dir, 'compartment/region4comp')
        system(sprintf('mkdir -p %s', region4comp_dir))

        #########################################################

        interSparse_ALL = foreach(chr_num2look=chr_nums, .combine=rbind) %dopar%
        {
            if(chr_num2look==23) chr_num2look = 'X'
            chr_name = paste0('chr', chr_num2look)
            chr2strawr = as.character(chr_num2look)

            if(is.null(green_list)) green_bins = as.character(1:1E7) ## all bins are 'green'
            if(!is.null(green_list))
            {  
                green_list_chr = green_list[chr==chr_name]
                if(nrow(green_list_chr) < (30E6 / bin_size)) return(NULL)
                green_bins = generate_binned_tab(input_GR=makeGRangesFromDataFrame(green_list_chr), bin_size)
                green_bins = as.character(green_bins$end / bin_size)
            }

            interDense_tab = foreach(mate=mates, .combine=merge) %do%
            {   
                hic_f = get_mate_hic(wk_dir, mate, impute=TRUE)             
                dump_mat = try(silent=TRUE, strawr::straw(matrix='observed', binsize=bin_size, chr1loc=chr2strawr, chr2loc=chr2strawr, unit='BP', fname=hic_f, norm='NONE'))
                if(class(dump_mat)[1]=="try-error") return(NULL)

                dump_mat = as.data.table(dump_mat)

                colnames(dump_mat) = c("pos_1", "pos_2", "counts")
                mat = convert_contact_format_fun(dump_mat, from='format_dump', to='format_sparse', chr_num=chr_num2look, binsize=bin_size)$format_sparse
                green_bins = intersect(rownames(mat), green_bins)
                mat = mat[green_bins, green_bins]

                g_mat = as.matrix(mat)
                g_mat[abs(col(mat) - row(mat)) > 10E6/bin_size] = NA
                n_max = nrow(g_mat)

                bins = 1:n_max
                inter_sparse = foreach(bin2look=1:n_max, .combine=c) %do%
                {
                    upper_indices = 1:bin2look
                    down_indieces = bin2look:n_max
                    inter_mat = g_mat[upper_indices, down_indieces]
                    inter_val = mean(inter_mat!=0, na.rm=TRUE) ## density
                    inter_val
                }

                bins = as.numeric(rownames(mat))

                tab = data.table(chr=chr_name, start=bins*bin_size, end=bins*bin_size + bin_size, val=inter_sparse)
                colnames(tab)[4] = mate
                tab
            }

            if(is.null(interDense_tab)) return(NULL)
            interDense_tab
        }

        interSparse_ALL = interSparse_ALL[mat > 0 & pat > 0]

        #########################################################

        interSparse_ALL$min_val = apply(interSparse_ALL[, mget(mates)], 1, min)
        interSparse_ALL$max_val = apply(interSparse_ALL[, mget(mates)], 1, max)

        chrs2keep = unique(interSparse_ALL$chr)
        interSparse_ALL = foreach(chr2look=chrs2keep, .combine=rbind) %do% ## normalize for possible copy number difference
        {
            tmp = interSparse_ALL[chr==chr2look]
            med = median(tmp$pat - tmp$mat)
            if(med < 0) tmp$pat = tmp$pat - med
            if(med > 0) tmp$mat = tmp$mat + med
            tmp
        }

        interSparse_ALL$diff = abs(log2(interSparse_ALL$pat / interSparse_ALL$mat))
        for(mate in c(mates, 'diff')) 
        {
            bed_f = sprintf('%s/inter_sparsity.%s.bedgraph', region4comp_dir, mate)
            write.table(interSparse_ALL[, mget(c(GR_3_cols, mate))], file=bed_f, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
            tab2write = interSparse_ALL[, mget(c(GR_3_cols, mate))]

            tab2write[[4]] = 1 - tab2write[[4]]

            bed_f = sprintf('%s/inter_sparsity4manuscript.%s.bedgraph', region4comp_dir, mate)
            write.table(tab2write, file=bed_f, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

        }

        #########################################################

        block2keep = foreach(chr_name=chrs2keep, .combine=rbind) %dopar%
        {
            interDense_tab = interSparse_ALL[chr==chr_name]
            thresh_diff2use = max(quantile(interDense_tab$diff, 0.75)*3, 1)

            # interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$mat = interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$pat = 0 ## if there is big difference between mat and pat, in the meantime, at least sparse in one condition
            interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$mat = interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$pat = 0 ## if there is big difference between mat and pat, in the meantime, at least sparse in one condition
            interDense_tab[diff > 1.3]$mat = interDense_tab[diff > 1.3]$pat = 0 ## if there is a really big difference between mat and pat

            interDense_tab$index = 1:nrow(interDense_tab)

            for(mate in mates) 
            {
                vec_mate = 1*(interDense_tab[[mate]] > thresh_mate)
                interDense_tab[[paste0('chunk_', mate)]] = c(0, cumsum(abs(diff(vec_mate))))
            }
        
            interDense_tab$chunk = paste0(interDense_tab$chunk_mat, interDense_tab$chunk_pat)
            chunk_li = split(interDense_tab, by='chunk')
            block2keep = chunk_li[[which.max(sapply(chunk_li, nrow))]][, mget(GR_3_cols)]
            if(nrow(block2keep) < 30E6 / bin_size) block2keep = NULL
            block2keep
        }

        block2keepGR = reduce(makeGRangesFromDataFrame(block2keep), min.gapwidth=5E6) ## remove small gaps, 16 Marh 2023 
        tab = as.data.table(block2keepGR)[, 1:3]
        tab$start = tab$start + 1
        block2keep = expand_comp_tab(tab, 'start', 'end', bin_size)[, c(1,5,6)]
        colnames(block2keep) = GR_3_cols

        bed_f = sprintf('%s/region4comp.bed', region4comp_dir)
        write.table(block2keep, file=bed_f, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        return(block2keep)
    }


    generate_comp_tab_pair = function(comp_tab_line, line_pair, chr_name=NA)
    {
        comp_pair = comp_tab_line[line_pair]
        if(!is.na(chr_name)) for(i in 1:length(comp_pair)) comp_pair[[i]] = comp_pair[[i]][chr==chr_name]
        for(i in 1:2) colnames(comp_pair[[line_pair[i]]])[3:5] = paste0(colnames(comp_pair[[line_pair[i]]])[3:5], "_", mates[i])
        comp_tab_pair = Reduce(merge, comp_pair)
        comp_tab_pair = comp_tab_pair[order(chr, bin_index)]

        comp_tab_pair = unique(comp_tab_pair, by=c("chr", "bin_index"))

        comp_tab_pair$comp_rank_diff = comp_tab_pair$comp_rank_pat - comp_tab_pair$comp_rank_mat
        comp_tab_pair$continous_rank_diff = comp_tab_pair$continous_rank_pat - comp_tab_pair$continous_rank_mat

        comp_tab_pair = get_flanking_diff(comp_tab_pair)
        comp_tab_pair$comp_name_mat = substr(comp_tab_pair$comp_name_mat, 1, 5)
        comp_tab_pair$comp_name_pat = substr(comp_tab_pair$comp_name_pat, 1, 5)
        comp_tab_pair$chr = as.character(comp_tab_pair$chr)
        return(comp_tab_pair)
    }

    #########################################################

    get_sig_comp_diff_thresh = function(comp_tab_line, lines_rand_combine, chr_num2look, comp_diff_id, quantile_threshs)
    {
        my_merge = function(...) merge(..., by=c("chr", "bin_index"))
        cols2choose = c("chr", "bin_index", comp_diff_id, "continous_rank_mat", "continous_rank_pat")
        chr_name2look = paste0('chr', chr_num2look)
        comp_diff_flank_rand = foreach(j=1:nrow(lines_rand_combine), .combine=my_merge) %dopar% 
        {
            comp_tab_pair_rand = generate_comp_tab_pair(comp_tab_line=comp_tab_line, unlist(lines_rand_combine[j,]), chr_name=chr_name2look)[, mget(cols2choose)]
            comp_tab_pair_rand$mean_comp = (comp_tab_pair_rand$continous_rank_mat + comp_tab_pair_rand$continous_rank_pat) / 2
            colnames(comp_tab_pair_rand)[c(3,6)] = c(sprintf("comp_diff_flank_%s", j), sprintf("mean_comp_%s", j))
            comp_tab_pair_rand[, c(3,6)] = abs(comp_tab_pair_rand[, c(3,6)])
            comp_tab_pair_rand[, c(1,2,3,6)]
        }

        diff_cols = colnames(comp_diff_flank_rand)[grepl("comp_diff_flank", colnames(comp_diff_flank_rand))]
        mean_comp_cols = colnames(comp_diff_flank_rand)[grepl("mean_comp", colnames(comp_diff_flank_rand))]

        sig_shift_thresh = data.table(sig_thresh=t(apply(comp_diff_flank_rand[, mget(diff_cols)], 1, function(v) quantile(v, quantile_threshs))))

        sig_shift_thresh$mean_comp = apply(comp_diff_flank_rand[, mget(mean_comp_cols)], 1, mean)
        sig_shift_thresh = cbind(comp_diff_flank_rand[, c("chr", "bin_index")], sig_shift_thresh)
        sig_thresh_cols = colnames(sig_shift_thresh)[grepl("sig_thresh", colnames(sig_shift_thresh))]


        sig_shift_thresh = sig_shift_thresh[, mget(c("chr", "bin_index", sig_thresh_cols))] 
        return(sig_shift_thresh)
    }


    get_flanking_diff = function(comp_tab_pair)
    {
        li = split(comp_tab_pair, by="chr")

        get_flanking_diff_helper = function(pat_vec, mat_vec)
        {
            bins2shift = 2
            shift_upstream = foreach(k=1:bins2shift, .combine=cbind) %do% 
            {
                diff_vec = pat_vec - c(mat_vec[-(1:k)], rep(0, k))
                data.table(diff_vec)
            }

            shift_no = foreach(k=0, .combine=cbind) %do% 
            {
                diff_vec = pat_vec - mat_vec ## no shift
                data.table(diff_vec)
            }

            shift_downstream = foreach(k=1:bins2shift, .combine=cbind) %do% 
            {
                n = length(mat_vec)
                diff_vec = pat_vec - c(rep(0, k), mat_vec[-((n-k+1):n)])
                data.table(diff_vec)
            }

            diff_shifts = cbind(shift_upstream, shift_no, shift_downstream)
            diff_flank = apply(diff_shifts, 1, function(v) {v[which.min(abs(v))]})
            return(diff_flank)  
        }

        n_chrs = length(li)
        for(i in 1:n_chrs) li[[i]]$comp_diff_flank = get_flanking_diff_helper(li[[i]]$comp_rank_pat, li[[i]]$comp_rank_mat)
        for(i in 1:n_chrs) li[[i]]$continous_rank_diff_flank = get_flanking_diff_helper(li[[i]]$continous_rank_pat, li[[i]]$continous_rank_mat)
        
        return(do.call( rbind, li ))
    }

    ##########################################################################################

    get_bad_comp_chr = function(LINE) ## chrs with failed compartment calls
    {
        load(unil_comp_lines_cor_file)
        lines_ALL_raw = expand.grid( c(LINE, paste0(LINE, '_random_k=', 1:5)), mates )
        lines_ALL = paste0(lines_ALL_raw[,1], '.', lines_ALL_raw[,2])
        bad_comp_chrs_raw = foreach(LINE=lines_ALL, .combine=c) %dopar%
        {
            cor_tab = comp_tab_lines_cor[[LINE]]
            cor_tab_ALL = merge(data.table(chr_num=1:22), cor_tab, all.x=TRUE)
            cor_tab_ALL$cor[is.na(cor_tab_ALL$cor)] = 0
            bad_comp_chr = cor_tab_ALL[cor < 0.4]$chr_num
            if(grepl('HCT_116', LINE) | grepl('HCT116', LINE) ) bad_comp_chr = cor_tab_ALL[cor < 0.1]$chr_num
            bad_comp_chr
        }

        if(LINE=='NHA') bad_comp_chrs_raw = setdiff(bad_comp_chrs_raw, 13)
        if(grepl('RPE', LINE) & !grepl('pWGD_T', LINE)) bad_comp_chrs_raw = setdiff(bad_comp_chrs_raw, 13) ## pWGD_T (Tumor): chr13 one copy lost
        return(sort(unique(bad_comp_chrs_raw)))
    }


    ################################ this script identifies low coverage region, and will be used to remove such regions

    get_sparse_and_low_coverage_region = function(wk_dir, bin_size)
    {
        # bin_size = 50E3 ## changed to 50kb, 6 Oct 2023 / remove any 40kb
        CNV_dir = file.path(wk_dir, sprintf("CNV/%skb/", bin_size/1E3))
        dir.create(CNV_dir, recursive=TRUE)
        for(map_type in c("phased", 'unphased')) compute_coverage_sparsity(wk_dir, map_type, bin_size=bin_size)
        chr_GR = makeGRangesFromDataFrame(generate_region_tab(bin_size=bin_size)) 

        ################################

        pdf_f = file.path(CNV_dir, sprintf("coverage.hist.pdf"))
        pdf(pdf_f, width=8, height=10)
        par(mfrow=c(2,1))

        low_coverage_GR = foreach(mate=mates, .combine=c) %do% 
        {
            coverage_file = file.path(CNV_dir, sprintf("coverage.%s.bedgraph", mate))
            coverage_tab = fread(coverage_file)[, 1:4]
            coverage_tab[[4]] = log2(coverage_tab[[4]] + 1)
            colnames(coverage_tab) = c(GR_3_cols, 'coverage')

            coverage_thresh = min(7, max(quantile(coverage_tab[coverage > 0]$coverage, 0.05), 6)) ## remove the lowest 5%. Set the thresh as 6 to be more conservative / reduce false positives
            # if(ID %in% IDs.RPE) coverage_thresh = max(quantile(coverage_tab[coverage > 0]$coverage, 0.05), 4.5) ## remove the lowest 10%

            main = sprintf("Contact matrix sparsity (or coverage) plot | thresh=%s", substr(coverage_thresh, 1, 4))
            hist(coverage_tab$coverage, breaks=100, xlim=c(0, 15), col=rgb(1,0,0,0.4), xlab="mat", ylab="", main=main); abline(v=coverage_thresh, col='black', lty=2, lwd=2)
            coverage_low = coverage_tab[coverage < coverage_thresh]
            low_GR = reduce(makeGRangesFromDataFrame(coverage_low), min.gapwidth=1E6)
            low_coverage_bed = file.path(CNV_dir, sprintf("low_coverage.%s.bedgraph", mate))
            fwrite(as.data.table(low_GR), file=low_coverage_bed, sep='\t', quote=FALSE, col.names=FALSE)
            low_GR
        }

        dev.off()

        low_coverage_GR = reduce(low_coverage_GR, min.gapwidth=1E6)
        chrs_blank = setdiff(chr_chars[1:23], as.character(seqnames(low_coverage_GR)))
        blank_chr_GR = chr_GR[seqnames(chr_GR) %in% chrs_blank]
        low_coverage_GR = c(low_coverage_GR, blank_chr_GR)

        low_coverage_bed = file.path(CNV_dir, "low_coverage.bed")
        fwrite(as.data.table(low_coverage_GR), file=low_coverage_bed, sep='\t', quote=FALSE, col.names=FALSE)

        ################################

        pdf_f = file.path(CNV_dir, sprintf("sparsity.hist.pdf"))
        pdf(pdf_f, width=8, height=10)
        par(mfrow=c(2,1))

        high_sparsity_GR = foreach(mate=mates, .combine=c) %do% 
        {
            sparsity_file = file.path(CNV_dir, sprintf("sparsity.%s.bedgraph", mate))
            sparsity_tab = fread(sparsity_file)[, 1:4]
            colnames(sparsity_tab) = c(GR_3_cols, 'sparsity')
            sparsity_tab$sparsity = 1 - sparsity_tab$sparsity

            # sparsity_thresh = min(0.9, quantile(sparsity_tab[sparsity < 1]$sparsity, 0.95)) ## remove the lowest 10%
            # if(ID %in% IDs.RPE) sparsity_thresh = min(0.95, quantile(sparsity_tab[sparsity < 1]$sparsity, 0.98)) ## remove the lowest 10%

            sparsity_thresh = 0.95 ## mainly use coverage to filter

            main = sprintf("Contact matrix sparsity (or sparsity) plot | thresh=%s", substr(sparsity_thresh, 1, 4))
            hist(sparsity_tab$sparsity, breaks=100, xlim=c(0.4, 1), col=rgb(1,0,0,0.4), xlab="mat", ylab="", main=main); abline(v=sparsity_thresh, col='black', lty=2, lwd=2)
            sparsity_high = sparsity_tab[sparsity > sparsity_thresh]
            high_GR = reduce(makeGRangesFromDataFrame(sparsity_high), min.gapwidth=1E6)
            high_sparsity_bed = file.path(CNV_dir, sprintf("high_sparsity.%s.bedgraph", mate))
            fwrite(as.data.table(high_GR), file=high_sparsity_bed, sep='\t', quote=FALSE, col.names=FALSE)
            high_GR
        }

        dev.off()

        high_sparsity_GR = reduce(high_sparsity_GR, min.gapwidth=1E6)
        high_sparsity_bed = file.path(CNV_dir, "high_sparsity.bed")
        fwrite(as.data.table(high_sparsity_GR), file=high_sparsity_bed, sep='\t', quote=FALSE, col.names=FALSE)
        

        ################################ combine the two / black list region

        black_list_region = reduce(c(high_sparsity_GR, low_coverage_GR), min.gapwidth=1E6)
        black_list_region_bed = file.path(CNV_dir, sprintf("black_list.bed"))

        ################################ setdiff to get green region

        ref = 'hg19'
        chr_size = fread(chrom_sizes)
        colnames(chr_size) = c("chr", "len")
        chr_size = chr_size[chr %in% paste0("chr", c(1:22, "X"))]
        chr_size_GR = makeGRangesFromDataFrame(data.table(chr=chr_size$chr, start=1, end=chr_size$len))
        end(chr_size_GR) = end(chr_size_GR) %/% bin_size * bin_size

        green_region_GR = GenomicRanges::setdiff(chr_size_GR, black_list_region)
        green_region_bed = file.path(CNV_dir, sprintf("green_list.bed"))

        ################################

        black_list_region = generate_binned_tab(black_list_region, bin_size=bin_size)
        green_region_tab = generate_binned_tab(green_region_GR, bin_size=bin_size)

        for(i in 2:3) black_list_region[[i]] = as.character(black_list_region[[i]]) 
        for(i in 2:3) green_region_tab[[i]] = as.character(green_region_tab[[i]]) 

        fwrite(black_list_region, file=black_list_region_bed, sep='\t', quote=FALSE, col.names=FALSE)
        fwrite(green_region_tab, file=green_region_bed, sep='\t', quote=FALSE, col.names=FALSE) 
    }


	diffComp_fun = function(wk_dir, chr_nums2look)
	{
		chr_names2look = paste0('chr', chr_nums2look)
		my_merge = function(...) merge(..., by=c("chr", "bin_index"))
		comp_diff_id = c("comp_diff_flank", "continous_rank_diff_flank")[2]
		bin_size = 10E3
		line_pair = sprintf("phased.%s", mates)

		lines = c('phased', sprintf("random_k=%s", 1:5)) ## do not include line without random
		comp_line_f = file.path(wk_dir, "/compartment/comp_tab_line.Rdata")
		comp_line_cor_f = file.path(wk_dir, "/compartment/comp_tab_line.cor.Rdata")

		################################

		load(comp_line_f); load(comp_line_cor_f)
		for(i in 1:length(comp_tab_line)) comp_tab_line[[i]][['chr']] = as.character(comp_tab_line[[i]][['chr']])	

		################################

		comp_tab_pair = generate_comp_tab_pair(comp_tab_line=comp_tab_line, line_pair=line_pair)
		colnames(comp_tab_pair)[colnames(comp_tab_pair)==comp_diff_id] = "comp_diff"

		################################ get_region with valid (random) compartment annotation

		comp_chr_bins = Reduce(my_merge, lapply(comp_tab_line, function(v) v[, 1:2]))[chr %in% chr_names2look]
		comp_chr_bins$start = (comp_chr_bins$bin_index - 1)*10E3 + 1
		comp_chr_bins$end = comp_chr_bins$bin_index*10E3
		GR = makeGRangesFromDataFrame(comp_chr_bins)
		
		bed_raw = as.data.table( reduce(GR, min.gapwidth=5E3) )[, 1:3]
		colnames(bed_raw) = c('chr', 'pos_start', 'pos_end')
		pos_start = as.character(bed_raw$pos_start)
		pos_end = as.character(bed_raw$pos_end)

		# comp_green_bed = data.table(chr=paste0('chr', bed_raw$chr), pos_start, pos_end, '.', 1, '.', pos_start, pos_end, '#f3c99d')
		# green_list_region_bed = file.path(unil_result_dir, sprintf("4.CNV/%s/%s.comp_green_list.bed", ID, ID)) ## reliable compartment region
		# try(fwrite(comp_green_bed, file=green_list_region_bed, sep='\t', quote=FALSE, col.names=FALSE))

		# ################################

		(lines_rand = names(comp_tab_line)[grepl("rand", names(comp_tab_line))])
		lines_rand_combine = data.table( t(combn(lines_rand, 2)) )
		
		#################################

		brut_wig_file = file.path(wk_dir, 'compartment/diffComp/diffComp.brut.wig')
		sig_wig_file = file.path(wk_dir, 'compartment/diffComp/diffComp.sig.wig')
		rand_wig_file = file.path(wk_dir, 'compartment/diffComp/diffComp.rand.wig')
		sig_seg_file = file.path(wk_dir, 'compartment/diffComp/diffComp.sig.bed')

		unlink(brut_wig_file)
		unlink(sig_wig_file)
		unlink(rand_wig_file)
		unlink(sig_seg_file)

		#################################

		quantile_thresh = "sig_thresh.80%"
		quantile_threshs = c(0.5, 0.8, 0.9, 0.95)

		sigDiffComp_bed_ALL = foreach(chr_num2look=chr_nums2look, .combine=rbind) %do%
		{
			cat(chr_num2look, sep="\n")
			sig_comp_change_helper(comp_tab_line, comp_tab_pair=comp_tab_pair, lines_rand_combine=lines_rand_combine, chr_num2look, quantile_threshs, quantile_thresh, comp_diff_id)
		}

		sigDiffComp_bed_ALL[sigDiffComp_bed_ALL=='chr23'] = 'chrX'
		fwrite(sigDiffComp_bed_ALL, file=sig_wig_file, sep="\t", col.names=FALSE)

		# sigDiffComp_bed_lines[[ID]] = sigDiffComp_bed_ALL
		# save(sigDiffComp_bed_lines, file=sigDiffComp_bed_lines_file) ## do not save, due to using %dopar%, Y.L, 02 Aug 2022

		################################# define consecutive segment
		
		directions = c("positive", "negative")

		#################################
		
		sigDiffComp_seg = foreach(direction=directions, .combine=rbind) %dopar%
		{
			seg_tab = foreach(name=unique(sigDiffComp_bed_ALL$chr), .combine=rbind) %do%
			{
				seg_tab_raw = sigDiffComp_bed_ALL[chr==name]
				if(direction=="negative") seg_tab_raw$score = -seg_tab_raw$score	

				seg_tab_raw = seg_tab_raw[score > 0.05] ## at bin level, only keep those > 0.05
				if(nrow(seg_tab_raw) == 0) return(NULL)

				seg_tab_raw$bin_index = as.numeric(seg_tab_raw$pos_end) / 10E3
				seg_GR = reduce(min.gapwidth=40E3, makeGRangesFromDataFrame(data.table(chr=seg_tab_raw$chr, start=seg_tab_raw$pos_start, end=seg_tab_raw$pos_end, score=seg_tab_raw$score)))

				seg_info_tab = foreach(k=1:length(seg_GR), .combine=rbind) %do%
				{
					bins = ceiling((start(seg_GR[k])/10E3)) : (end(seg_GR[k]) / 10E3)
					mean_seg_score = mean(seg_tab_raw[bin_index %in% bins]$score)
					data.table(n_bins=length(bins), mean_score=mean_seg_score)
				}

				seg_tab = data.table(as.data.table(seg_GR)[, 1:3], seg_info_tab)
				colnames(seg_tab)[1:3] = c("chr", "pos_start", "pos_end")
				seg_tab = seg_tab[n_bins > 10 & mean_score > 0.1] ## longer than 10 bins & at seg level, only keep those > 0.1
				seg_tab = seg_tab[n_bins > 30 | mean_score > 0.2] ## add more strengent threshold on 27 Feb 2022, selected based on the correlation and consistency of HCT116 with H3K9me3
				seg_tab
			}

			for(i in 2:3) seg_tab[[i]] = as.character(seg_tab[[i]])
			sigDiffComp_seg = data.table(seg_tab[, 1:3], substr(seg_tab$mean_score, 1, 4), seg_tab$mean_score, ".", seg_tab[, 2:3], col=ifelse(direction=="positive", "#FF0000", "#0000FF"))
			sigDiffComp_seg
		}

		#################################

		fwrite(sigDiffComp_seg, file=sig_seg_file, sep="\t", col.names=FALSE)	
	}

	sig_comp_change_helper = function(comp_tab_line, comp_tab_pair, lines_rand_combine, chr_num2look, quantile_threshs, quantile_thresh, comp_diff_id)
	{
		###################################################### enhancer peak

		chr_name2look = paste0("chr", chr_num2look)
		sig_comp_shift_thresh = get_sig_comp_diff_thresh(comp_tab_line, lines_rand_combine, chr_num=chr_num2look, comp_diff_id=comp_diff_id, quantile_threshs=quantile_threshs)

		comp_diff_pair = merge(comp_tab_pair[chr==chr_name2look], sig_comp_shift_thresh, by=c("chr", "bin_index"))
		comp_diff_pair$is_comp_diff_sig = 1*(abs(comp_diff_pair$comp_diff) > comp_diff_pair[[quantile_thresh]])
		comp_diff_pair$diff = sign(comp_diff_pair$comp_diff) * (abs(comp_diff_pair$comp_diff) - 1*comp_diff_pair[[quantile_thresh]])

		# brut_comp_shift_wig = data.table(chr=chr_name2look, pos_start=as.character(comp_diff_pair$pos_start), pos_end=as.character(comp_diff_pair$pos_end), score=comp_diff_pair$comp_diff) ## the brut comp diff, without deduction
		# rand_comp_shift_wig = data.table(chr=chr_name2look, pos_start=as.character((sig_comp_shift_thresh$bin_index - 1)*10E3), pos_end=as.character(sig_comp_shift_thresh$bin_index*10E3), score=sig_comp_shift_thresh[[quantile_thresh]])

		# fwrite(brut_comp_shift_wig, file=brut_wig_file, sep="\t", col.names=FALSE, append=TRUE)
		# fwrite(rand_comp_shift_wig, file=rand_wig_file, sep="\t", col.names=FALSE, append=TRUE)

		sigDiffComp = comp_diff_pair[is_comp_diff_sig==1]
		pos_start = as.character(sigDiffComp$pos_start)
		pos_end = as.character(sigDiffComp$pos_end)
		sigDiffComp_bed = data.table(chr=sigDiffComp$chr, pos_start=pos_start, pos_end=pos_end, score=sigDiffComp$diff)	
		return(sigDiffComp_bed)
	}
