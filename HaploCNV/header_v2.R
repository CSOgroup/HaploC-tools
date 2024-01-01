
    
    ## This file generates the input for diploid.pl
    ## Yuanlong LIU

    #########################################################

    get_encode_black_bins = function(bin_size) ## black list region from encode
    {
        black_fs = file.path(HaploCNV_dir, '/ENCODE_blacklist/', c('ENCFF001TDO.bed', 'ENCFF001THR.bed', 'hg19-blacklist.v2.bed'))
        black_tab = fread(black_fs[3])[, 1:3] ## the third one seems most relevant
        colnames(black_tab) = GR_3_cols
        black_GR = makeGRangesFromDataFrame(black_tab)
        black_GR = black_GR + 0E3
        black_bins = generate_binned_tab(black_GR, bin_size=bin_size)
        return(black_bins)
    }


    get_lines2use = function(return_tab=FALSE)
    {
        post_status_f = file.path(post_script_dir, '7.Post_analysis/11.Post_status_lines/post_status_lines.tsv')
        post_status_tab = fread(post_status_f)
        post_status_tab = post_status_tab[!(ID %in% c('GM19099', 'KBM7', 'GSE161493_LNCaP'))] ## GM19099: low phased hic; KBM7: 1 ploidy; GSE161493_LNCaP: not profiled
        post_status_tab = post_status_tab[!grepl('RIF1AID', ID, ignore.case=TRUE)]
        LINEs = post_status_tab$ID
        cat('n_LINEs:', length(LINEs), '\n')
        if(return_tab) return(post_status_tab)
        return(LINEs)
    }

    convert_ID2Name2show = function(IDs)
    {
        post_status_tab = get_lines2use(return_tab=TRUE)
        Name2show = post_status_tab$Name2show[ match(IDs, post_status_tab$ID) ]
        return(Name2show)
    }    
    
    #########################################################


    chr_nums = c(1:22, 'X')
    chr_chars = paste0('chr', c(1:22, 'X'))

    chr_pairs_intra = data.table(expand.grid(1:23, 1:23))[Var1 == Var2]
    chr_pairs_intra[chr_pairs_intra=='23'] = 'X'

    chr_pairs_inter = data.table(expand.grid(1:23, 1:23))[Var1 < Var2]
    chr_pairs_inter[chr_pairs_inter=='23'] = 'X'

    ######################################################### retrieve biggest block to compute compartments

    get_black_region = function(LINE) ## region with sparse values should be removed from both mat and pat
    {
        black_list_region_f = sprintf("%s/4.CNV/%s/%s.high_sparsity.bed", unil_result_dir, LINE, LINE)
        black_list_region = fread(black_list_region_f)[, 1:3]
        colnames(black_list_region) = GR_3_cols
        black_list_region = as.data.table(reduce(makeGRangesFromDataFrame(black_list_region), min.gapwidth=1E6))[, 1:3]
        colnames(black_list_region) = GR_3_cols
        black_list_region$chr = as.vector(black_list_region$chr)
        return(black_list_region)
    }


    # get_largest_block = function(bin_size=100E3, LINE, green_list=NULL, chr_chars, thresh_mate) ## if density < 0.2 / largest block for deriving compartment / no much need to correct for CNV, because the criteria seems powerful to handle CNV / 25-02-2023
    # {
    #     region4comp_dir = file.path(unil_calder_dir, LINE, 'region4comp')
    #     system(sprintf('mkdir -p %s', region4comp_dir))

    #     chr_chars = setdiff(chr_chars, 'chrX')
    #     chr_nums = as.numeric(gsub('chr', '', chr_chars))

    #     interSparse_ALL = foreach(chr_num_x=chr_nums, .combine=rbind) %dopar%
    #     {
    #         chr_name = paste0("chr", chr_num_x)   
    #         chr2strawr = as.character(chr_num_x)
    #         if(is.null(green_list)) green_bins = as.character(1:1E7) ## all bins are 'green'
    #         if(!is.null(green_list))
    #         {  
    #             green_list_chr = green_list[chr==chr_name]
    #             if(nrow(green_list_chr) < (30E6 / bin_size)) return(NULL)
    #             green_bins = generate_binned_tab(input_GR=makeGRangesFromDataFrame(green_list_chr), bin_size)
    #             green_bins = as.character(green_bins$end / bin_size)
    #         }

            
    #         interDense_tab = foreach(mate=mates, .combine=merge) %do%
    #         {
    #             mat_dump_deleted = strawr::straw(matrix='observed', binsize=bin_size, chr1loc=chr2strawr, chr2loc=chr2strawr, unit='BP', fname=get_mate_hic(LINE, mate), norm='NONE') %>% data.table
    #             colnames(mat_dump_deleted) = c("pos_1", "pos_2", "counts")
    #             mat = convert_contact_format_fun(mat_dump_deleted, from='format_dump', to='format_sparse', chr_num=chr_num_x, binsize=bin_size)$format_sparse
    #             green_bins = intersect(rownames(mat), green_bins)
    #             mat = mat[green_bins, green_bins]

    #             g_mat = as.matrix(mat)
    #             g_mat[abs(col(mat) - row(mat)) > 10E6/bin_size] = NA
    #             n_max = nrow(g_mat)

    #             bins = 1:n_max
    #             inter_sparse = foreach(bin2look=1:n_max, .combine=c) %do%
    #             {
    #                 upper_indices = 1:bin2look
    #                 down_indieces = bin2look:n_max
    #                 inter_mat = g_mat[upper_indices, down_indieces]
    #                 inter_val = mean(inter_mat!=0, na.rm=TRUE) ## density
    #                 inter_val
    #             }

    #             bins = as.numeric(rownames(mat))

    #             tab = data.table(chr=chr_name, start=bins*bin_size, end=bins*bin_size + bin_size, val=inter_sparse)
    #             colnames(tab)[4] = mate
    #             tab
    #         }

    #         if(is.null(interDense_tab)) return(NULL)
    #         interDense_tab
    #     }

    #     interSparse_ALL = interSparse_ALL[mat > 0 & pat > 0]

    #     #########################################################

    #     interSparse_ALL$min_val = apply(interSparse_ALL[, mget(mates)], 1, min)
    #     interSparse_ALL$max_val = apply(interSparse_ALL[, mget(mates)], 1, max)

    #     chrs2keep = unique(interSparse_ALL$chr)
    #     interSparse_ALL = foreach(chr2look=chrs2keep, .combine=rbind) %do% ## normalize for possible copy number difference
    #     {
    #         tmp = interSparse_ALL[chr==chr2look]
    #         med = median(tmp$pat - tmp$mat)
    #         if(med < 0) tmp$pat = tmp$pat - med
    #         if(med > 0) tmp$mat = tmp$mat + med
    #         tmp
    #     }

    #     interSparse_ALL$diff = abs(log2(interSparse_ALL$pat / interSparse_ALL$mat))
    #     for(mate in c(mates, 'diff')) 
    #     {
    #         bed_f = sprintf('%s/inter_sparsity.%s.bedgraph', region4comp_dir, mate)
    #         write.table(interSparse_ALL[, mget(c(GR_3_cols, mate))], file=bed_f, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    #         tab2write = interSparse_ALL[, mget(c(GR_3_cols, mate))]

    #         tab2write[[4]] = 1 - tab2write[[4]]

    #         bed_f = sprintf('%s/inter_sparsity4manuscript.%s.bedgraph', region4comp_dir, mate)
    #         write.table(tab2write, file=bed_f, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

    #     }

    #     #########################################################

    #     block2keep = foreach(chr_char=chrs2keep, .combine=rbind) %dopar%
    #     {
    #         interDense_tab = interSparse_ALL[chr==chr_char]
    #         thresh_diff2use = max(quantile(interDense_tab$diff, 0.75)*3, 1)

    #         # interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$mat = interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$pat = 0 ## if there is big difference between mat and pat, in the meantime, at least sparse in one condition
    #         interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$mat = interDense_tab[diff > thresh_diff2use & min_val < thresh_mate]$pat = 0 ## if there is big difference between mat and pat, in the meantime, at least sparse in one condition
    #         interDense_tab[diff > 1.3]$mat = interDense_tab[diff > 1.3]$pat = 0 ## if there is a really big difference between mat and pat

    #         interDense_tab$index = 1:nrow(interDense_tab)

    #         for(mate in mates) 
    #         {
    #             vec_mate = 1*(interDense_tab[[mate]] > thresh_mate)
    #             interDense_tab[[paste0('chunk_', mate)]] = c(0, cumsum(abs(diff(vec_mate))))
    #         }
        
    #         interDense_tab$chunk = paste0(interDense_tab$chunk_mat, interDense_tab$chunk_pat)
    #         chunk_li = split(interDense_tab, by='chunk')
    #         block2keep = chunk_li[[which.max(sapply(chunk_li, nrow))]][, mget(GR_3_cols)]
    #         if(nrow(block2keep) < 30E6 / bin_size) block2keep = NULL
    #         block2keep
    #     }

    #     block2keepGR = reduce(makeGRangesFromDataFrame(block2keep), min.gapwidth=5E6) ## remove small gaps, 16 Marh 2023 
    #     tab = as.data.table(block2keepGR)[, 1:3]
    #     tab$start = tab$start + 1
    #     block2keep = expand_comp_tab(tab, 'start', 'end', bin_size)[, c(1,5,6)]
    #     colnames(block2keep) = GR_3_cols

    #     bed_f = sprintf('%s/region4comp.bed', region4comp_dir)
    #     write.table(block2keep, file=bed_f, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    #     return(block2keep)
    # }
    

    ######################################################### gc_content, 03 Nov 2022

    get_GC = function(bin_size=100E3) ## http://bioconductor.org/help/course-materials/2014/SeattleOct2014/A01.2_IntroductionToBioconductor.html
    {
        suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

        chr_size = foreach(chr_char=chr_chars, .combine=rbind) %dopar% data.table(chr=chr_char, len=length(BSgenome.Hsapiens.UCSC.hg19[[chr_char]]))
        chr_size = chr_size[chr %in% paste0("chr", c(1:22, "X"))]

        region_tab = foreach(chr_query=chr_size$chr, .combine=rbind) %dopar%
        {
            n_bins = chr_size[chr==chr_query]$len %/% bin_size - 1
            region_tab_chr = data.table(chr=chr_query, start=(0:n_bins)*bin_size + 1)
            region_tab_chr$end = region_tab_chr$start + bin_size - 1
            region_tab_chr
        }

        region_GR = makeGRangesFromDataFrame(region_tab)
        seqsInRegion <- getSeq(BSgenome.Hsapiens.UCSC.hg19, region_GR)
        alf1 <- as.data.table(alphabetFrequency(seqsInRegion, collapse=FALSE))
        alf1$atgc = apply(alf1[, c('A', 'T', 'G', 'C')], 1, sum)
        alf1$gc = apply(alf1[, c( 'G', 'C')], 1, sum)        
        alf1$gc_ratio = alf1$gc / alf1$atgc
        gc_tab = cbind(region_tab, alf1[, 'gc_ratio'])
        return(gc_tab)
    }

    ######################################################### get CNV neutral region

    split_GR2bin = function(GR, bin_size) ## convert granges to bin
    {
        GR_li = split(GR, seqnames(GR))
        ref_binned = foreach(gr=GR_li, .combine=c) %do%
        {
            start = floor(min(start(gr)) / bin_size) * bin_size + 1
            end = ceiling(max(end(gr)) / bin_size) * bin_size
            gr_bins = slidingWindows(GRanges(seqnames(gr)[1], IRanges(start=min(start(gr)), end=max(end(gr)))), width=bin_size, step=bin_size)[[1]] 
            gr_bins
        }

        GR_binned = ref_binned[unique(queryHits(findOverlaps(ref_binned, GR)))]
        return(GR_binned)
    }

    get_cnv_bamlance_region = function(LINE, ref='hg19')
    {
        chr_size = fread(chrom_sizes)
        colnames(chr_size) = c("chr", "len")
        chr_size = chr_size[chr %in% paste0("chr", c(1:22, "X"))]
        chr_size_GR = makeGRangesFromDataFrame(data.table(chr=chr_size$chr, start=1, end=chr_size$len))

        ## cnv_segmentation
        segment_bed_f = sprintf("%s/4.CNV/%s/%s.my_segment.bedgraph", unil_result_dir, LINE, LINE)
        segments = fread(segment_bed_f)[, 1:3]
        colnames(segments) = GR_3_cols
        segments_GR = makeGRangesFromDataFrame(segments)
        segments_GR = reduce(segments_GR, min.gapwidth=100E3) ## if there is a gap of 100E3, these 100E3 might be difficult to analyze, exclude
        bamlance_GR = setdiff(chr_size_GR, segments_GR, ignore.strand=TRUE)
        return(bamlance_GR)
    }

    ######################################################### differential reads

    diffReads_fun = function(input_bam_mates, output_bed, region_tab, log2_thresh, n_max) ## count number of reads mapped to each bin
    {
        region_tab = region_tab[chr %in% paste0('chr', 1:22)]
        my_combine = function(...) merge(..., by=c('chr', 'pos_start', 'pos_end'))
        n_reads = foreach(mate=mates, .combine=my_combine) %do%
        {
            n_reads = featureCount(input_bam_mates[mate], output_bed=NULL, region_tab)
            colnames(n_reads)[4] = mate
            return(n_reads)
        }

        n_reads$log2FC = log2((n_reads$pat + 1) / (n_reads$mat + 1))
        n_reads$n_max = apply(n_reads[, mget(mates)], 1, max)
        n_max_thresh = n_max
        n_reads$pass_nmax = 1*(n_reads$n_max >= n_max_thresh)
        n_reads$pass_log2thresh = 1*(abs(n_reads$log2FC) > log2_thresh)  
        n_reads$n_max_thresh = n_max_thresh
        n_reads_ALL = n_reads[n_max >= n_max_thresh] ## testable regiosn  
        n_reads = n_reads[n_max >= n_max_thresh & abs(log2FC) > log2_thresh]
        fwrite(n_reads[, c(1:3, 6)], file=output_bed, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
        fwrite(n_reads_ALL[, c(1:3, 6)], file=gsub('log2.bedgraph', 'log2.ALL.bedgraph', output_bed), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
        return(n_reads_ALL)
    }

    #########################################################

 

    get_ins = function(LINE, mate='mat')
    {
        ins_f = sprintf("%s/3.Structural/2.TADs/%s/%s.ins.%s.bedgraph", unil_result_dir, LINE, LINE, mate)
        return(ins_f)
    }

    get_diffIns = function(LINE)
    {
        diffIns_f = sprintf("%s/3.Structural/2.TADs/%s/diffTAD.%s.thresh=0.01.bed", unil_result_dir, LINE, LINE)
        return(diffIns_f)
    } 

    get_diffComp = function(LINE)
    {
        diffComp_f = sprintf("%s/3.Structural/1.diffComp/%s/%s.diffComp.sig.bed", unil_result_dir, LINE, LINE)
        return(diffComp_f)
    } 

    get_opt_cor = function(LINE, mate)
    {
        opt_cor_f = sprintf("%s/Calder/%s/%s/sub_compartments/cor_with_ref.ALL.txt", wk_dir, LINE, mate)
        opt_cor_val_f = sprintf("%s/Calder/%s/%s/sub_compartments/cor_with_ref.txt", wk_dir, LINE, mate)
        opt_cor_fs = c(opt_cor_f=opt_cor_f, opt_cor_val_f=opt_cor_val_f)
        return(opt_cor_fs)
    }  

    get_rand_opt_cor = function(LINE, mate)
    {
        LINE = paste0(LINE, '_random_k=5')
        opt_cor_f = sprintf("%s/Calder/%s/%s/sub_compartments/cor_with_ref.ALL.txt", wk_dir, LINE, mate)
        opt_cor_val_f = sprintf("%s/Calder/%s/%s/sub_compartments/cor_with_ref.txt", wk_dir, LINE, mate)
        opt_cor_fs = c(opt_cor_f=opt_cor_f, opt_cor_val_f=opt_cor_val_f)
        return(opt_cor_fs)
    }    

    ######################################################### get_abs_cnv

    get_phased_cnv = function(totalCNV, pat_ratio) 
    {
        A <- matrix(c(1, 1- pat_ratio, 1, -pat_ratio), nrow=2)
        b <- c(totalCNV, 0)
        cnv_mates = as.vector(solve(A) %*% t(t(b)))
        names(cnv_mates) = c('mat', 'pat')
        return(cnv_mates)
    }

    #########################################################

    get_hap_config_tab = function()
    {
        ################################ common ratios

        hap_config_tab = as.data.table(expand.grid(n_mat=1:5, n_pat=1:5))
        hap_config_tab$n_total = hap_config_tab$n_mat + hap_config_tab$n_pat        
        hap_config_tab$n_mated = sprintf('%s|%s', hap_config_tab$n_mat, hap_config_tab$n_pat)
        hap_config_tab$ratio = hap_config_tab$n_pat / (hap_config_tab$n_mat + hap_config_tab$n_pat)
        hap_config_tab = hap_config_tab[n_pat >= n_mat]
        hap_config_tab = hap_config_tab[n_mat + n_pat < 6]
        hap_config_tab = unique(hap_config_tab, by='ratio')

        hap_config_tab$ratio_min = 0
        hap_config_tab$ratio_max = 1

        hap_config_tab[n_mated=='1|1']$ratio_min = 0.48 ## 0.500
        hap_config_tab[n_mated=='1|1']$ratio_max = 0.55

        hap_config_tab[n_mated=='1|2']$ratio_min = 0.62 ## 0.667
        hap_config_tab[n_mated=='1|2']$ratio_max = 0.67

        hap_config_tab[n_mated=='1|3']$ratio_min = 0.70 ## 0.750
        hap_config_tab[n_mated=='1|3']$ratio_max = 0.80

        hap_config_tab[n_mated=='1|4']$ratio_min = 0.75 ## 0.800
        hap_config_tab[n_mated=='1|4']$ratio_max = 0.85

        hap_config_tab[n_mated=='2|3']$ratio_min = 0.55 ## 0.600
        hap_config_tab[n_mated=='2|3']$ratio_max = 0.65

        hap_config_tab = hap_config_tab[order(n_mat, n_pat)]
        hap_config_tab$index = 1:nrow(hap_config_tab)

        hap_config_tab$col = comp_cols[c(1, 5:8)]
        return(hap_config_tab)
    }

    get_hap_config_tab_relaxed = function()
    {
        ################################ common ratios

        hap_config_tab = as.data.table(expand.grid(n_mat=1:15, n_pat=1:15))
        hap_config_tab$n_total = hap_config_tab$n_mat + hap_config_tab$n_pat        
        hap_config_tab$n_mated = sprintf('%s|%s', hap_config_tab$n_mat, hap_config_tab$n_pat)
        hap_config_tab$ratio = hap_config_tab$n_pat / (hap_config_tab$n_mat + hap_config_tab$n_pat)
        hap_config_tab = hap_config_tab[n_pat >= n_mat]
        hap_config_tab = unique(hap_config_tab, by='ratio')

        # hap_config_tab$ratio_min = 0
        # hap_config_tab$ratio_max = 1
        # hap_config_tab[n_mated=='0|x']$ratio_min = 0.85 ## 0.600
        # hap_config_tab[n_mated=='0|x']$ratio_max = 1

        hap_config_tab = hap_config_tab[order(n_mat, n_pat)]
        hap_config_tab$index = 1:nrow(hap_config_tab)
        hap_config_tab = hap_config_tab[ratio < 0.85] ## more likely to be x|0

        hap_config_tab$col = rep(comp_cols[c(1, 5:8)], length=nrow(hap_config_tab))
        return(hap_config_tab)
    }

    make_bed = function(tab, color_name='red', bed_f)
    {
        pos_start = as.character(tab$start)
        pos_end = as.character(tab$end)
        bed = data.table(tab$chr, pos_start, pos_end, 'loss', '1', '.', '1', pos_start, pos_end, col2hex(color_name))
        fwrite(bed, file=bed_f, sep='\t', quote=FALSE, col.names=FALSE)
    }

    expand_comp_tab = function(tab, start_col='V2', end_col='V3', bin_size=10E3)
    {
        bin_indices_li = apply(tab[, mget(c(start_col, end_col))], 1, function(v) {bins = 1 + seq(v[1]-1, v[2], by=bin_size)/bin_size; bins[1:(length(bins) - 1)]})
        if(class(bin_indices_li)[1]=='matrix') bin_indices_li = lapply(seq_len(ncol(bin_indices_li)), function(i) bin_indices_li[,i]) ## if equal width of each row

        expand_len = sapply(bin_indices_li, length)

        tab = tab[rep(1:nrow(tab), expand_len), ]
        tab$bin_index = unlist(bin_indices_li)
        tab$pos_start = (tab$bin_index - 1)*bin_size + 1
        tab$pos_end = (tab$bin_index - 0)*bin_size
        return(tab)
    }

    expand_comp_tab_v2 = function(tab, start_col='start', end_col='end', bin_size=10E3)
    {
        bin_indices_li = apply(tab[, mget(c(start_col, end_col))], 1, function(v) {bins = 1 + seq(v[1]-1, v[2], by=bin_size)/bin_size; bins[1:(length(bins) - 1)]})
        if(class(bin_indices_li)[1]=='matrix') bin_indices_li = lapply(seq_len(ncol(bin_indices_li)), function(i) bin_indices_li[,i]) ## if equal width of each row

        expand_len = sapply(bin_indices_li, length)

        tab = tab[rep(1:nrow(tab), expand_len), ]
        tab$bin_index = unlist(bin_indices_li)
        tab$start = (tab$bin_index - 1)*bin_size + 1
        tab$end = (tab$bin_index - 0)*bin_size
        return(tab)
    }

    ##############################

    get_hic_coverage_withInter = function(LINE, bin_size=40E3, mate, flanking_bp=NA, type='hic', impute=TRUE)
    {
        if(mate=='bulk') hic_f = get_bulk_hic(LINE)

        coverage_chrs_tmp = foreach(chr_num_x=1:23, .combine=rbind) %do%
        {
            cat(chr_num_x ,'\n')
            chr_char_x = paste0('chr', chr_num_x)
            contacts_chrs = foreach(chr_num_y=1:23, .combine=merge_all) %do%
            {
                mat_dump_deleted = try(silent=TRUE, strawr::straw(matrix='observed', binsize=bin_size, chr1loc=as.character(ifelse(chr_num_x==23, 'X', chr_num_x)), chr2loc=as.character(ifelse(chr_num_y==23, 'X', chr_num_y)), unit='BP', fname=hic_f, norm='NONE') %>% data.table)
                mat_dump_deleted = na.omit(mat_dump_deleted)

                if(chr_num_y < chr_num_x) colnames(mat_dump_deleted)[1:2] = colnames(mat_dump_deleted)[2:1]
                {
                    mat_dump_deleted[, sum:=sum(counts), by='x']
                    tmp = unique(mat_dump_deleted[order(x)], by='x')[, c('x', 'sum')]
                    tmp$chr_bin = paste0(chr_char_x, '_', tmp$x / bin_size + 1)
                    tmp = tmp[, 3:2]
                    colnames(tmp)[2] = sprintf('n_chr%s', chr_num_y)
                }
                tmp
            }

            contacts_chrs$contacts = apply(contacts_chrs[, -1], 1, sum, na.rm=TRUE)
            contacts_chrs[, c('chr_bin', 'contacts')]
        }

        ######################################
        
        chr = do.call(rbind, strsplit(coverage_chrs_tmp$chr, '_'))[,1]
        chr[chr=='chr23'] = 'chrX'
        pos_start = (as.numeric(do.call(rbind, strsplit(coverage_chrs_tmp$chr, '_'))[,2]) - 1)*bin_size + 1
        pos_end = pos_start + bin_size - 1

        coverage_chrs = data.table(chr=chr, start=pos_start, end=pos_end, contacts=coverage_chrs_tmp$contacts)
        coverage_chrs = coverage_chrs[order(chr, start)]
        for(z in 2:3) coverage_chrs[[z]] = as.character(coverage_chrs[[z]])

        bedgraph_file = file.path(unil_result_dir, sprintf("4.CNV/%s/%s.coverage.%s.inter.bedgraph", LINE, LINE, mate))

        fwrite(coverage_chrs, file=bedgraph_file, sep='\t', quote=FALSE, col.names=FALSE)
        return(coverage_chrs)
    }


        featureCount = function(input_bam, output_bed, region_tab, chr_chars=NULL) ## count number of reads mapped to each bin
        {
            if(is.null(chr_chars[1])) chr_chars = paste0('chr', 1:22)
            region_tab = region_tab[chr %in% chr_chars]
            gtf_header = "/mnt/etemp/Yuanlong/3.Packages/RSEM/my_tutorial_and_test/test/ref_hg19/gencode.v19.annotation.gtf"
            key = paste(sample(letters, 20), sep="", collapse="")
            gtf_tmp = sprintf("/mnt/etemp/Yuanlong/4.Tmp/delete_%s.gtf", key)
            feature_out_tmp = sprintf("/mnt/etemp/Yuanlong/4.Tmp/delete_%s.fc", key)

            region_tab[[2]] = as.character(region_tab[[2]])
            region_tab[[3]] = as.character(region_tab[[3]])


            options(scipen=20)

            gtf_tab = data.table(region_tab$chr, "enhancer", "exon", region_tab[, 2:3], ".", "+", ".", "gene_id \"")
            gtf_tab[[9]] = paste0(gtf_tab[[9]], "peak_", 1:nrow(gtf_tab), "\";")
            gtf_tab[nrow(gtf_tab), "V3"] = "exon"

            system(sprintf("head -n 5 %s > %s", gtf_header, gtf_tmp))
            fwrite(gtf_tab, file=gtf_tmp, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
            cmd = sprintf("~/miniconda3/envs/intg/bin/featureCounts -O -a %s -o %s %s -T 20 -g gene_id &> /dev/null", gtf_tmp, feature_out_tmp, input_bam) ## http://subread.sourceforge.net/featureCounts.html
            system(cmd)
            tab = fread(feature_out_tmp)
            bed = data.table(chr=tab$Chr, pos_start=as.character(tab$Start), pos_end=as.character(tab$End), count=tab[[7]])
            if(!is.null(output_bed)) fwrite(bed, file=output_bed, row.names=FALSE, col.names=FALSE, sep="\t")
            unlink(feature_out_tmp)
            unlink(gtf_tmp)
            return(bed)
        }


        generate_region_tab = function(bin_size, ref='hg19')
        {
            chr_size = fread(chrom_sizes)
            colnames(chr_size) = c("chr", "len")
            chr_size = chr_size[chr %in% paste0("chr", c(1:22, "X"))]

            region_tab = foreach(chr_query=chr_size$chr, .combine=rbind) %dopar%
            {
                n_bins = chr_size[chr==chr_query]$len %/% bin_size
                region_tab_chr = data.table(chr=chr_query, pos_start=(0:n_bins)*bin_size + 1)
                region_tab_chr$pos_end = region_tab_chr$pos_start + bin_size - 1
                region_tab_chr
            }
            return(region_tab)
        }

        generate_binned_tab = function(input_GR, bin_size)
        {
            chr_GR = makeGRangesFromDataFrame(generate_region_tab(bin_size))
            overlap_tab = as.data.table(chr_GR[unique(queryHits(findOverlaps(chr_GR, input_GR)))])[, 1:3]
            colnames(overlap_tab) = GR_3_cols
            overlap_tab$chr = as.vector(overlap_tab$chr)
            return(overlap_tab)
        }

        ####################################### generate chr_pos and mat_pat for shapeit
       

        generate_gene_coordinates = function(bin_size)
        {
            try(suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene)))
            try(suppressMessages(library(org.Hs.eg.db)))
            gene_info = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

            gene_symbol_raw <- AnnotationDbi::select(org.Hs.eg.db, keys=gene_info$gene_id, columns="SYMBOL", keytype="ENTREZID")
            gene_symbol = gene_symbol_raw$SYMBOL
            names(gene_symbol) = gene_symbol_raw$ENTREZID

            gene_cordinates = data.frame(chr=as.vector(seqnames(gene_info)), gene_id=gene_info$gene_id, pos_start=as.numeric(start(gene_info)), pos_end=as.numeric(end(gene_info)), gene_len=as.numeric(width(gene_info)), strand=strand(gene_info))
                
            pstrand_indices = which(gene_cordinates$strand == '+')
            nstrand_indices = which(gene_cordinates$strand == '-')
            gene_cordinates$TSS = 0
            gene_cordinates$TSS[pstrand_indices] = gene_cordinates[pstrand_indices, 'pos_start']
            gene_cordinates$TSS[nstrand_indices] = gene_cordinates[nstrand_indices, 'pos_end']
            gene_cordinates$TSS_bin_index = ceiling(gene_cordinates$TSS / bin_size)
            gene_cordinates$SYMBOL = unname(gene_symbol[gene_cordinates$gene_id])
            gene_cordinates = data.table(subset(gene_cordinates, chr %in% paste0("chr", c(1:22, "X"))))
            gene_cordinates = gene_cordinates[!is.na(SYMBOL)]
            return(gene_cordinates)
        }

        get_dump_mat = function(save_dir, chr_num, norm, bin_size, ob_oe=c('ob', 'oe')[1])  ## this will return a symetric mat
        {
            identity = "bulk"
            chr_num_x = chr_num_y = chr_num
            out_file = sprintf('%s/mat_chr%s_chr%s_%skb_%s.%s.%s.txt.gz', save_dir, chr_num_x, chr_num_y, bin_size/1E3, norm, ob_oe, identity)
            dump_tab = fread(out_file)
            colnames(dump_tab) = c("pos_1", "pos_2", "counts")
            dump_tab = dump_tab[!is.na(counts) & !is.nan(counts)]
            contact_mat = convert_contact_format_fun(input_dat=dump_tab, from='format_dump', to="format_sparse", chr_num=chr_num, binsize=bin_size)$format_sparse
            return(contact_mat)
        }

    generate_comp_tab_pair = function(line_pair, chr_num=NA)
    {
        comp_pair = comp_tab_line[line_pair]
        if(!is.na(chr_num)) for(i in 1:length(comp_pair)) comp_pair[[i]] = comp_pair[[i]][chr==chr_num]
        for(i in 1:2) colnames(comp_pair[[line_pair[i]]])[3:5] = paste0(colnames(comp_pair[[line_pair[i]]])[3:5], "_", mates[i])
        comp_tab_pair = Reduce(merge, comp_pair)
        comp_tab_pair = comp_tab_pair[order(chr, bin_index)]
        # comp_tab_pair$bin_index = ceiling( comp_tab_pair$bin_index / 2 ) ## 10kb to 20kb / if only I am studying hicdc

        comp_tab_pair = unique(comp_tab_pair, by=c("chr", "bin_index"))

        comp_tab_pair$comp_rank_diff = comp_tab_pair$comp_rank_pat - comp_tab_pair$comp_rank_mat
        comp_tab_pair$continous_rank_diff = comp_tab_pair$continous_rank_pat - comp_tab_pair$continous_rank_mat

        comp_tab_pair = get_flanking_diff(comp_tab_pair)
        comp_tab_pair$comp_name_mat = substr(comp_tab_pair$comp_name_mat, 1, 5)
        comp_tab_pair$comp_name_pat = substr(comp_tab_pair$comp_name_pat, 1, 5)
        return(comp_tab_pair)
    }


   

    plot_colorByDensity = function(x1,x2,
                               ylim=c(min(x2),max(x2)),
                               xlim=c(min(x1),max(x1)),
                               xlab="", ylab="", main="", cex=2) {
     
    df <- data.frame(x1,x2)
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
    df$col <- cols[df$dens]
    plot(x2~x1, data=df[order(df$dens),], 
         ylim=ylim,xlim=xlim,pch=20,col=col,
         cex=cex,xlab=xlab,ylab=ylab,
         main=main)
    }



    densplot2 = function(x, y, points = FALSE, pch=19, cex=1, xlim=c(min(x), max(x)), ylim=c(min(y), max(y)), reg_line=FALSE, x_thresh=0, y_thresh=0, loess_span=0.75, lm_fun=NULL, ...){
            df = data.table(x=x, y=y)
            df = df[order(x)]
            x = df$x
            y = df$y
            d = densCols(x, y, colramp=colorRampPalette(c("black", "white")))
            df$dens = col2rgb(d)[1,] + 1L

            cols = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)

            df$col = cols[df$dens]
            # df=df[order(dens)]
            if(points)
            points(df$x, df$y, pch=pch, col=df$col, cex=cex, ...)
            else
            plot(df$x, df$y, pch=pch, col=df$col, cex=cex, xlim=xlim, ylim=ylim, ...)
            
            if(lm_fun=='lm') 
            {
                abline(lm(df$y ~ df$x), col='blue')
                return(0)
            }


            if(reg_line)
            {
                df = unique(df[abs(x) >= x_thresh & abs(y) >= y_thresh])
                ls_line <- lm_fun(df$y ~ df$x)

                if(lm_fun==lm) plot()

                plx <- predict(ls_line, se=T)

                lines(ls_line$x, ls_line$fitted, col="red", lwd=2)
                lines(ls_line$x, ls_line$fitted - qt(0.975, plx$df)*plx$se, lty=2)
                lines(ls_line$x, ls_line$fitted + qt(0.975, plx$df)*plx$se, lty=2)

                x_long = c(ls_line$x, rev(ls_line$x))
                y_long = c(ls_line$fitted - qt(0.975, plx$df)*plx$se, rev(ls_line$fitted + qt(0.975, plx$df)*plx$se))
                polygon(x_long, y_long, col=rgb(0,0,1, 0.2), border=NA)
            }
    }



    get_black_list_bins = function(black_list_tab, bin_size) ## black_list_tab: chr, pos_start, pos_end; bin_size: the bin_size to be converted to
    {
        colnames(black_list_tab) = c('chr', 'pos_start', 'pos_end')
        black_list_tab_chrs = split(black_list_tab, by ='chr')
        black_bins_chrs = foreach(chr_name=names(black_list_tab_chrs)) %do%
        {
            black_list_chr_GR = makeGRangesFromDataFrame(black_list_tab_chrs[[chr_name]])
            max_pos = max(GenomicRanges::end(black_list_chr_GR)) + bin_size
            pos_start = seq(1, max_pos, by=bin_size)
            pos_end = pos_start + bin_size - 1

            bins_GR = makeGRangesFromDataFrame(data.table(chr=chr_name, pos_start=pos_start, pos_end=pos_end))

            overlaps = findOverlaps(bins_GR, black_list_chr_GR, minoverlap=10) ## 10 is some arbitary number
            black_bins = as.character( unique(queryHits(overlaps)) )
            return(black_bins)
        }
        names(black_bins_chrs) = names(black_list_tab_chrs)
        return(black_bins_chrs)
    }


    exclude_index_in_region = function(query_region, LINE)
    {
        black_list_f = sprintf('%s/4.CNV/%s/%s.black_list.bed', unil_result_dir, LINE, LINE)
        black_region = fread(black_list_f)[, 1:3]
        colnames(black_region) = GR_3_cols
        black_region_GR = reduce(makeGRangesFromDataFrame(black_region), min.gapwidth=100E3)
        query_region_GR = makeGRangesFromDataFrame(query_region)
        indices2rm = unique(queryHits(findOverlaps(query_region_GR, black_region_GR)))
        return(indices2rm)
    }
    
    get_contact_mat_fun = function(ratio2keep, chr_nums) ## keep_ratio = - down_ratio ## inter-chr
    {
        chr_names = paste0('chr', chr_nums)
        output_raw = foreach(i=chr_nums) %dopar%
        {
            inner_output = foreach(j=chr_nums) %dopar%
            {
                mat_file = paste0(mat_dir, 'mat_chr', j, '_chr', i, '_50kb_ob.txt')
                if(!file.exists(mat_file)) mat_file = paste0(mat_dir, 'mat_chr', i, '_chr', j, '_50kb_ob.txt')
                mat = fread(mat_file)

                if(ratio2keep!=1)
                {
                    counts_total = sum(mat[[3]])
                    counts2keep = ceiling(counts_total*ratio2keep)
                    set.seed(1)
                    samples = sample(1:nrow(mat), counts2keep, replace=TRUE, prob=mat[[3]])
                    freq = table(samples)
                    # freq = freq[order(as.numeric(names(freq)))]
                    mat = mat[as.numeric(names(freq)),]
                    mat[[3]] = as.vector(unname(freq))
                    mat = mat[order(V1, V2)]
                }

                mat[[1]] = mat[[1]] / resolution + 1
                mat[[2]] = mat[[2]] / resolution + 1

                oe_size_row = max(mat[[1]])
                oe_size_col = max(mat[[2]])

                mat_sparse = Matrix(0, nrow=oe_size_row, ncol=oe_size_col)
                row_col_size = max(oe_size_col, oe_size_row)

                if(i==j) mat_sparse = Matrix(0, nrow=row_col_size, ncol=row_col_size)
                mat_sparse[cbind(mat[[1]], mat[[2]])] = mat[[3]]
                if(i==j) mat_sparse = my_forceSymmetric(mat_sparse, uplo='U')

                rownames(mat_sparse) = 1:nrow(mat_sparse)
                colnames(mat_sparse) = 1:ncol(mat_sparse)
                mat_sparse
            }
            names(inner_output) = chr_names
            inner_output
        }

        #######
        names(output_raw) = chr_names
        output = output_raw
        get_chr_sizes = function(output)
        {
            chr_sizes = sapply(chr_names, function(chr) max(dim(output_raw[[chr]][[chr]])))
            return(chr_sizes)
        }

        chr_sizes = get_chr_sizes(output)

        for(i in chr_names)
        for(j in chr_names)
        {
            n_row = nrow(output_raw[[i]][[i]])
            n_col = ncol(output_raw[[i]][[i]])
            min_index = which.min(abs(chr_sizes[i] - c(n_row, n_col)))
            if(min_index==2 & (i<j)) output[[i]][[j]] = t(output[[i]][[j]])
            if(min_index==1 & (i>j)) output[[i]][[j]] = t(output[[i]][[j]])
        }

        if(length(chr_names)>1) chr_sizes_common = apply(sapply(chr_names, function(i) sapply(output[[i]], ncol)), 1, min)
        if(length(chr_names)==1) chr_sizes_common = ncol(output[[1]][[1]])
        names(chr_sizes_common) = chr_names

        for(i in chr_names)
        for(j in chr_names)
        {
            output[[i]][[j]] = output[[i]][[j]][1:chr_sizes_common[i], 1:chr_sizes_common[j]]
            rownames(output[[i]][[j]]) = paste0(i, '_', 1:chr_sizes_common[i])
            colnames(output[[i]][[j]]) = paste0(j, '_', 1:chr_sizes_common[j])
        }
        res = do.call(rbind, lapply( chr_names, function(i) as.matrix(do.call(cbind, output[[i]]) )))
        return(res)
    }


    vcf2tab = function(vcf_file, gz=FALSE, phased=TRUE, save2wig=FALSE, consider_ps=FALSE, QUAL_thresh=-1, wig_ground_truth=NA, format_shapeit2_hap=FALSE, block_ratio_thresh=NA, which_blocks=1, vcf_tab_only=TRUE)
    {
        library(data.table)

        if(format_shapeit2_hap==TRUE)
        {
            wig_file = gsub('.haps', '.wig', vcf_file, ignore.case=TRUE)
            my_skip = 0
            vcf_tab = data.table::fread(vcf_file, sep=' ', header=FALSE)
            vcf_tab = vcf_tab[, c(1,3:7)]
            colnames(vcf_tab)[c(1,2)] = c('chr', 'POS')         
            vcf_tab$hap = paste0(vcf_tab[[5]], '|', vcf_tab[[6]]) 
            wig_diff_file = gsub('.haps', '.diff.wig', vcf_file, ignore.case=TRUE)
            n_phased = nrow(vcf_tab)
            vcf_tab = vcf_tab[, mget(c('chr', 'POS', 'hap'))]
        }

        if(!format_shapeit2_hap)
        {
            if(gz==TRUE)
            {
                wig_file = gsub('.vcf.gz', '.wig', vcf_file, ignore.case=TRUE)
                my_skip = as.numeric(strsplit(system(intern=TRUE, sprintf('zgrep -n "#CHROM\tPOS" %s', vcf_file)), ':')[[1]][1]) - 1
            }

            if(gz==FALSE)
            {
                wig_file = gsub('[.]vcf', '.wig', vcf_file, ignore.case=TRUE)
                my_skip = as.numeric(strsplit(system(intern=TRUE, sprintf('grep -n "#CHROM\tPOS" %s', vcf_file)), ':')[[1]][1]) - 1
            }

            wig_diff_file = gsub('.vcf', '.diff.wig', vcf_file, ignore.case=TRUE)

            vcf_tab = data.table::fread(vcf_file, sep='\t', skip=my_skip)
            vcf_tab = vcf_tab[, mget(setdiff(colnames(vcf_tab), 'INFO'))]
            vcf_tab = vcf_tab[QUAL > QUAL_thresh]
        }

        sample_name = tail(colnames(vcf_tab), 1)
        vcf_tab$hap = substr(vcf_tab[[sample_name]], 1, 3)
        vcf_tab = vcf_tab[hap %in% c('0/1', '0|1', '1|0')]

        ## If wig_ground_truth is not NA
        if(length(wig_ground_truth) > 1) vcf_tab = vcf_tab[POS %in% wig_ground_truth$pos_end]

        colnames(vcf_tab)[c(1,2)] = c('chr', 'POS') 

        vcf_tab_input = vcf_tab
        n_total_input_SNPs =  nrow(vcf_tab_input) ## total number of variants in the VCF file input for phasing

        #########################################################

        if(consider_ps) 
        {
            vcf_tab$PS = unname(sapply(vcf_tab[[sample_name]], function(v) tail(strsplit(v, ':')[[1]],1)))
            vcf_tab = vcf_tab[PS!='.']

            vcf_tab[, `:=` (block_ratio = .N/nrow(vcf_tab)), by = PS]
            vcf_tab = vcf_tab[order(-block_ratio, POS)]
            vcf_tab$block_index = c(1, cumsum(diff(vcf_tab$block_ratio)!=0) + 1)
            vcf_tab[, `:=` (block_range = paste0(substr(min(POS)/1E6, 1, 4), '_', substr(max(POS)/1E6, 1, 4))), by = PS]

            n_phased = nrow(vcf_tab)
            max_ps = names(which.max(table(vcf_tab$PS)))
            max_2block_total_coverage = sum(tail(sort(table(vcf_tab$PS)),2)) / n_phased

            if(is.na(block_ratio_thresh)) block_ratio_thresh = max(vcf_tab$block_ratio)
            if(is.na(which_blocks[1])) vcf_tab = vcf_tab[block_ratio >= block_ratio_thresh]
            if(!is.na(which_blocks[1]) & length(which_blocks) > 1) vcf_tab = vcf_tab[block_index %in% which_blocks]

            max_block_total_coverage = nrow(vcf_tab[block_ratio==max(vcf_tab$block_ratio)]) / n_total_input_SNPs
            max_block_phased_coverage = nrow(vcf_tab[block_ratio==max(vcf_tab$block_ratio)]) / n_phased
            selected_block_phased_coverage = nrow(vcf_tab) / n_phased

            cat('  >> n_max_block in n_total_input_SNPs:', max_block_total_coverage, ' | n_max_block in n_phased_SNPs:', max_block_phased_coverage, '| selected_blocks in n_phased_SNPs:', selected_block_phased_coverage, '\n')

        }

        hap_val = c(1,-1, 1)
        names(hap_val) = c('0|1', '1|0', '0/1')
        vcf_tab$hap_val = unname(hap_val[vcf_tab$hap])/(vcf_tab$block_index)

        if(save2wig==FALSE) return(vcf_tab)

        if(grepl('chr', vcf_tab[[1]])[1]) wig_tab = data.table(chr=vcf_tab[[1]], pos_start=vcf_tab$POS - 1, pos_end=vcf_tab$POS, hap=vcf_tab$hap_val)
        if(!grepl('chr', vcf_tab[[1]])[1]) wig_tab = data.table(chr=paste0('chr', vcf_tab[[1]]), pos_start=vcf_tab$POS - 1, pos_end=vcf_tab$POS, hap=vcf_tab$hap_val)

        # wig_tab$hap = hap_val[wig_tab$hap]
        wig_tab = na.omit(wig_tab)

        data.table::fwrite(wig_tab, file=wig_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

        if(!phased) return(vcf_tab)

        phased_ground_truth = consist_ratio = NA
        if(length(wig_ground_truth) > 1)
        {
            colnames(wig_ground_truth)[4] = 'hap_ground_truth'
            wig_overlap = merge(wig_tab, wig_ground_truth, by=c('chr', 'pos_start', 'pos_end'))

            phased_ground_truth = nrow(wig_overlap) / nrow(wig_ground_truth)

            snp_coverage = sum(vcf_tab_input$POS %in% wig_ground_truth$pos_end) / nrow(wig_ground_truth)
            phased_coverage = nrow(wig_overlap) / nrow(wig_ground_truth)

            consist_ratio_raw = mean(wig_overlap$hap==wig_overlap$hap_ground_truth)
            consist_ratio = max(consist_ratio_raw, 1 - consist_ratio_raw)
            cat('  >> ground truth snp coverage:', snp_coverage, '\n')
            cat('  >> ground truth phased coverage:', phased_coverage, ' | consist_ratio:', consist_ratio, ' | ')
            consist_sign = sign(consist_ratio_raw - 0.5)
            wig_overlap$hap = consist_sign*(wig_overlap$hap)
            wig_tab_diff = wig_overlap[hap!=hap_ground_truth]
            data.table::fwrite(wig_tab_diff, file=wig_diff_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
        }

        cat('done\n')
        if(consider_ps==TRUE) phasing_info = list(phased_ground_truth=phased_ground_truth, phased_in_total=n_phased/n_total_input_SNPs, max_block_ratio=max_block_phased_coverage, max_2block_ratio=max_2block_total_coverage, consist_ratio=consist_ratio, vcf_tab=vcf_tab)
        if(consider_ps==FALSE) phasing_info = list(phased_ground_truth=phased_ground_truth, max_block_ratio=NA, max_2block_ratio=max_2block_ratio, consist_ratio=consist_ratio, vcf_tab=vcf_tab)
        if(vcf_tab_only) phasing_info = vcf_tab
        return(phasing_info)
    }
    