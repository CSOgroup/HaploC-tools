
    trim_outlier = function(vec) ## 2*exp(2.5): copy > 24
    {
        vec[vec > 2.5] = 2.5; vec[vec < -2.5] = -2.5; return(vec)
    }

    my_mean = function(vec) ## do not consider the lowest 5%
    {
        # vec = vec[vec < (median(vec))]
        vec = vec[vec > 0]
        vec = vec[between(vec, quantile(vec, 0.05), quantile(vec, 1))]
        return(median(vec))
    }

    # get_deletions = function(ID, bin_size=100E3)
    # {
    #     bedgraph_f = file.path(unil_result_dir, sprintf("4.CNV/%s/%s.sparsity.%s.bedgraph", ID, ID, 'bulk'))
    #     sparsity_tab = fread(bedgraph_f)[, 1:4]
    #     colnames(sparsity_tab) = c(GR_3_cols, 'sparsity')
    #     del_region = sparsity_tab[sparsity > 0.65]
    #     if(bin_size==20E3) del_region = sparsity_tab[sparsity > 0.9]

    #     bed_f = sprintf('%s/%s.del_region.bed', mark_dir, ID)
    #     fwrite(del_region[, 1:3], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
    #     return(del_region[, 1:3])
    # }

    get_deletions = function(wk_dir) ## 10E3
    {
        mark_dir = file.path(wk_dir, 'CNV/10kb')
        bedgraph_f = file.path(mark_dir, sprintf("coverage.bulk.bedgraph"))
        sparsity_tab = fread(bedgraph_f)[, 1:4]
        colnames(sparsity_tab) = c(GR_3_cols, 'sparsity')
   
        problem_region = get_encode_black_bins(bin_size=10E3) ## problematic region to exclude
        sparsity_tab = sparsity_tab[!problem_region, on=GR_3_cols] 

        del_region = sparsity_tab[sparsity < 200]
        bed_f = sprintf('%s/del_region.bed', mark_dir)
        fwrite(del_region[, 1:3], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
        return(del_region[, 1:3])
    }

    add_comp = function(wk_dir) ## only region have compartment will be included, can update to include all regions
    {
        get_comp = function()
        {
            comp_tab_f = file.path(wk_dir, 'compartment/comp_tab_line.Rdata')
            load(comp_tab_f)
            comp = comp_tab_line[['bulk']][, c('chr', 'pos_start', 'pos_end', 'continous_rank')]
            colnames(comp) = c(GR_3_cols, 'comp_rank')
            # comp$chr = paste0('chr', comp$chr)
            comp = comp[end %% bin_size2look==0]
            comp$start = comp$end - bin_size2look + 1
            return(comp)
        } 
        comp = get_comp()
        coverage_mates = merge(coverage_mates, comp)
        coverage_mates = coverage_mates[order(chr, start)]
        return(coverage_mates)
    }

    ################################

    add_tech_features = function(coverage_mates, bin_size=100E3, extream=0.05)
    {
        rescale_fun = function(vec) ## trim extream values
        {
            vec[vec < quantile(vec, extream)] = quantile(vec, extream)
            vec[vec > quantile(vec, 1 - extream )] = quantile(vec, 1 - extream)
            return(vec)
        }
        bin_size_char = sprintf('%skb', bin_size/1E3)
        cutSitesNumber_gc_mapscore_tab = readRDS(cutSitesNumber_gc_mapscore_rds_f)
        cutSitesNumber_gc_mapscore_tab = cutSitesNumber_gc_mapscore_tab[[bin_size_char]]
        cutSitesNumber_gc_mapscore_tab[, bin_size:=NULL]

        cutSitesNumber_gc_mapscore_tab

        coverage_mates = merge(coverage_mates, cutSitesNumber_gc_mapscore_tab)

        for(key in c('gcper', 'mapscore', 'cutSitesNumber')) coverage_mates[[key]] = rescale_fun(coverage_mates[[key]])
        return(coverage_mates)
    }

    ################################ round to integer CNV based on whether the region was phased or not

    round_to_nearest <- function(x) ifelse(x %% 1 <= 0.5, floor(x), ceiling(x))

    # my_dev_fun = function(x, y)
    # {
    #     mean_xy = (x + y) /2
    #     out = (max(x, y) - mean_xy) / mean_xy
    #     return(out)
    # }


    my_dev_fun = function(x, y)
    {
        out = abs(log(x) - log(y))
        return(out)
    }

    ################################ find regions with low unphased / phased SNPs

    generate_snp_count_per_bin = function(wk_dir, bin_size=100E3)
    {
        snp_count_tab = foreach(chr2look=chr_chars2look, .combine=rbind) %do%
        {
            snp_tab_both = foreach(case=c('unphased', 'phased'), .combine=merge_all) %do%
            {
                if(case=='unphased')
                {
                    vcf_raw = sprintf('%s/mega/VCFs/mega.%s.HET.bi.QUAL=10.vcf', wk_dir, chr2look)
                    snp_tab = vcf2tab(vcf_raw, save2wig=FALSE)[, c('chr', 'POS')]              
                }

                if(case=='phased')
                {
                    vcf_phased = sprintf('%s/mega/HapCut/phased_hap_intg/hic_%s_opt.pat_mat', wk_dir, chr2look)
                    if(!file.exists(vcf_phased)) snp_tab = data.table(chr=chr2look, POS=1)
                    if(file.exists(vcf_phased))
                    {
                        snp_tab = fread(vcf_phased, header=FALSE)
                        snp_tab$chr = chr2look
                        snp_tab$POS = as.numeric(do.call(rbind, strsplit(snp_tab[[1]], ':'))[,2])                         
                    }
                }

                snp_tab = snp_tab[, c('chr', 'POS')]
                snp_tab$bin_index = ceiling(snp_tab$POS / bin_size)

                snp_tab[, n_snp_by_bin:=length(POS), by=bin_index]
                snp_tab = unique(snp_tab, by='bin_index')
                colnames(snp_tab)[4] = case
                snp_tab[, POS:=NULL]
                return(snp_tab)
            }

            snp_tab_both[is.na(snp_tab_both)] = 0

            start = (snp_tab_both$bin_index -1)*bin_size + 1
            end = snp_tab_both$bin_index*bin_size
            count_tab = data.table(chr=chr2look, start=start, end=end, n_unphased=snp_tab_both$unphased, n_phased=snp_tab_both$phased)
            count_tab
        }

        tab2write = snp_count_tab[, c('chr', 'start', 'end', 'n_unphased')]
        snp_count_f1 = sprintf('%s/mega/VCFs/unphased_snp_density.bedgraph', wk_dir)
        write.table(tab2write, file=snp_count_f1, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

        snp_count_f2 = sprintf('%s/CNV/unphased_snp_density.bedgraph', wk_dir)
        write.table(tab2write, file=snp_count_f2, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

        tab2write = snp_count_tab[, c('chr', 'start', 'end', 'n_phased')]
        snp_count_f2 = sprintf('%s/CNV/phased_snp_density.bedgraph', wk_dir)
        write.table(tab2write, file=snp_count_f2, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

        return(snp_count_tab)
    }

    get_low_snp_region = function(wk_dir, redo_count=FALSE) ## find regions with low SNP density. These regions will not be used for defining 1|1 and 1|2 CRCs
    {
        if(redo_count) generate_snp_count_per_bin(wk_dir, bin_size=100E3)
        chr_bins = generate_region_tab(100E3)
        colnames(chr_bins) = GR_3_cols

        snp_count_f = sprintf('%s/CNV/phased_snp_density.bedgraph', wk_dir)
        snp_low_f = sprintf('%s/CNV/low_snp_region.bed', wk_dir)

        snp_count_tab = fread(snp_count_f)
        colnames(snp_count_tab) = c(GR_3_cols, 'n_snp')
        low_thresh = 20
        high_region = snp_count_tab[n_snp > low_thresh]

        low_snp_region = chr_bins[!high_region[, 1:3], on = GR_3_cols]
        low_GR = makeGRangesFromDataFrame(low_snp_region)
        low_GR_merge = reduce(low_GR, min.gapwidth=0.5E6)
        low_GR_merge = low_GR_merge[width(low_GR_merge) > 2E6] ## more than 2Mb
        low_GR_merge$low_snp_seg = sprintf('%s:%s-%sMb', seqnames(low_GR_merge), (start(low_GR_merge) - 1)/1E6, end(low_GR_merge)/1E6)
        low_tab = generate_binned_tab(low_GR_merge, bin_size=100E3)
        low_snp_seg = low_GR_merge$low_snp_seg[subjectHits(findOverlaps(makeGRangesFromDataFrame(low_tab), low_GR_merge))]
        write.table(low_tab, file=snp_low_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
        low_tab$low_snp_seg = low_snp_seg
        return(low_tab)
    }

    ################################ define diffCNV_seg. 

    generate_diffCNV_seg = function(wk_dir, bin_size, region2exclude, low_snp_region)
    {
        region2exclude$end = (region2exclude$end %/% 100E3 + 1) * 100E3 ## 100E3: fixed size
        region2exclude$start = region2exclude$end - 100E3 + 1
        region2exclude = unique(region2exclude)

        low_snp_region$end = (low_snp_region$end %/% 100E3 + 1) * 100E3
        low_snp_region$start = low_snp_region$end - 100E3 + 1
        low_snp_region = unique(low_snp_region)        

        ################################

        mark_dir = file.path(wk_dir, "CNV")

        ################################

        coverage_pat_f = sprintf("%s/100kb/coverage.ratio.pat.bedgraph", mark_dir)
        segment_bed_f = sprintf("%s/%skb/diffCNV.bedgraph", mark_dir, bin_size/1E3)

        coverage_raw = fread(coverage_pat_f) ## work on paternal coverage
        colnames(coverage_raw) = c(GR_3_cols, 'ratio')
        coverage_raw$percent = coverage_raw$ratio + 0.5
        col2keep = colnames(coverage_raw)
        coverage_raw = coverage_raw[!region2exclude, on = GR_3_cols] 

        ################################

        coverage_new = foreach(chr=c(1:22, 'X'), .combine=rbind) %do% ## 
        {
            chr_char = paste0('chr', chr)
            coverage = coverage_raw[chr==chr_char][order(start)]
            if(nrow(coverage) < 100) return(NULL)

            coverage_n_low = coverage[!low_snp_region, on = GR_3_cols] ## not low snp region
            # coverage_low =  coverage[low_snp_region, on = GR_3_cols]  ## low snp region. There can still have meaningful diffCNV segmentation because of imputation
            coverage_low =  merge(coverage, low_snp_region)  ## low snp region

            coverage_n_low$index = 1:nrow(coverage_n_low)
            coverage_low$index = 1:nrow(coverage_low)

            cnv_tree <- rpart(ratio ~ index, data=coverage_n_low, control=rpart.control(minsplit=5, minbucket=5, cp=0.01)) ## segmentation
            # cnv_tree <- rpart(ratio ~ index, data=coverage_n_low, control=rpart.control(minsplit=1, minbucket=1, cp=0.01)) ## segmentation

            coverage_n_low$seg_val = predict(cnv_tree, data.table(index=coverage_n_low$index))
            coverage_n_low[, seg_var:=mean(abs(ratio - seg_val)), by='seg_val'] ## variation inside each segment

            coverage_low[, seg_val:=mean(ratio), by='low_snp_seg']
            coverage_low[, seg_var:=mean(abs(ratio - seg_val)), by='low_snp_seg'] ## variation inside each segment
            coverage_low[, low_snp_seg:=NULL]

            coverage = rbind(coverage_n_low, coverage_low)[order(start)]
            # coverage = coverage[seg_var < 0.05 & abs(seg_val) < ((4/(4+1)) - 0.5)] ## exclude high var regions
            coverage = coverage[seg_var < 0.1 & abs(seg_val) < ((4/(4+1)) - 0.5)] ## exclude high var regions
            return(coverage[, mget(col2keep)])
        }
        
        coverage_seg_tab = foreach(chr=c(1:22, 'X'), .combine=rbind) %do% ## 
        {
            chr_char = paste0('chr', chr)
            coverage = coverage_new[chr==chr_char][order(start)]
            if(nrow(coverage) < 100) return(NULL)

            coverage$index = 1:nrow(coverage)

            cnv_tree <- rpart(ratio ~ index, data=coverage, control=rpart.control(minsplit=5, minbucket=5, cp=0.02)) ## segmentation
            coverage$seg_val = predict(cnv_tree, data.table(index=coverage$index))
            coverage[, seg_var:=mean(abs(ratio - seg_val)), by='seg_val'] ## variation inside each segment

            coverage = coverage[seg_var < 0.05 & abs(seg_val) < ((4/(4+1)) - 0.5)] ## exclude high var regions
            return(coverage)
        }

        coverage_seg_tab = expand_comp_tab_v2(coverage_seg_tab, bin_size=bin_size)
        
        make_seg_bed = function(tab, bed_f)
        {
            start = as.character(tab$start)
            end = as.character(tab$end)
            bed = data.table(tab$chr, start, end, tab$seg_val)
            fwrite(bed, file=bed_f, sep='\t', quote=FALSE, col.names=FALSE)
        }

        make_seg_bed(coverage_seg_tab, segment_bed_f)
    }

    ################################ low snp region setdiff generate_diffCNV_seg region are considered as LOH regions 

    add_LOH_info = function(coverage_mates_final, wk_dir, bin_size)
    {
        LOH_candidate = coverage_mates_final[rCN=='n/p'][, mget(GR_3_cols)] ## regions without proper diffCNV, these can be LOH or unphased region
        snp_count_f1 = sprintf('%s/mega/VCFs/unphased_snp_density.bedgraph', wk_dir)
        snp_count_tab = fread(snp_count_f1)
        colnames(snp_count_tab) = c(GR_3_cols, 'n_snp')

        snp_count_tab = expand_comp_tab_v2(snp_count_tab, bin_size=bin_size)
        snp_count_tab[, bin_index:=NULL]

        LOH_candidate_GR = makeGRangesFromDataFrame(LOH_candidate)
        GRs = reduce(LOH_candidate_GR)
        LOH_candidate$GR_index = subjectHits(findOverlaps(LOH_candidate_GR, GRs))

        LOH_candidate = merge(LOH_candidate, snp_count_tab)
        LOH_candidate[, median_snp:=median(1.0*n_snp), by='GR_index']
        LOH = LOH_candidate[median_snp < 30]

        LOH_f = sprintf('%s/CNV/LOH_region.bed', wk_dir)
        write.table(LOH[, mget(GR_3_cols)], file=LOH_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

        LOH$LOH = 1
        coverage_mates_final = merge(coverage_mates_final, LOH, all.x=TRUE)
        coverage_mates_final[is.na(LOH)]$LOH = 0
        return(coverage_mates_final)
    }

   set_aCN_11_12 = function(coverage_mates, which2use, version)
    {
        cat('set absolute CN in 1:1 1:2 region | ')
        coverage_mates$which2use = coverage_mates[[which2use]]
        CRC_11_12 = foreach(categ=c('1|1', '1|2'), .combine=rbind) %do%
        {
            coverage_mates_cat = coverage_mates[rCN==categ]
            coverage_mates_cat[, cov_by_chunk:=1.0*median(which2use), by='chunk']
            coverage_mates_cat[, chunk_size:=length(bulk)*bin_size2look/1E6, by='chunk']
            uniq_chunk_tab = unique(coverage_mates_cat[, c('chunk', 'cov_by_chunk', 'chunk_size')])
            uniq_chunk_tab = uniq_chunk_tab[chunk_size > 10] ## consider chunk with size > 10Mb ## set this value large
            uniq_chunk_tab$rCN = categ
            uniq_chunk_tab
        }

        CRC_11_12$tCN = 1
        CRC_11_12$cov_by_chunk = CRC_11_12$cov_by_chunk / min(CRC_11_12$cov_by_chunk) / 1.1

        is_11_min = CRC_11_12$rCN[which.min(CRC_11_12$cov_by_chunk)]=='1|1'
        is_12_min = !is_11_min

        if(is_11_min)
        {
            cat('bulk min comes from: 11\n')
            base = 2
            CRC_11_12$cov_by_chunk = base*CRC_11_12$cov_by_chunk
            for(i in 1:nrow(CRC_11_12))
            {
                if(CRC_11_12$rCN[i]=='1|1') candidates = c(2, 4, 6, 8,  10)
                if(CRC_11_12$rCN[i]=='1|2') candidates = c(3, 6, 9, 12, 15)
                tCN = candidates[which.min(abs(CRC_11_12$cov_by_chunk[i] - candidates))]
                CRC_11_12$tCN[i] = tCN
            }
        }

        if(is_12_min)
        {
            cat('Warning: bulk min comes from: 12\n')

            # bulk_min_from_cat_f = sprintf('%s/4.CNV/bulk_min_from_cat.Rds', unil_result_dir)
            # bulk_min_from_cat_li = readRDS(bulk_min_from_cat_f)
            # bulk_min_from_cat_li[[ID]] = '1|2'
            # saveRDS(bulk_min_from_cat_li, file=bulk_min_from_cat_f)

            base = 3
            CRC_11_12$cov_by_chunk = base*CRC_11_12$cov_by_chunk
            for(i in 1:nrow(CRC_11_12))
            {
                if(CRC_11_12$rCN[i]=='1|1') candidates = c(4, 6, 8, 10)
                if(CRC_11_12$rCN[i]=='1|2') candidates = c(3, 6, 9, 12)
                tCN = candidates[which.min(abs(CRC_11_12$cov_by_chunk[i] - candidates))]
                CRC_11_12$tCN[i] = tCN
            }
        }

        for(j in 1:nrow(CRC_11_12)) coverage_mates[chunk==CRC_11_12$chunk[j]]$tCN = CRC_11_12$tCN[j] ## update with new tCN
        
        ################################

        coverage_mates[, which2use:=NULL]
        bed_f = sprintf('%s/tCN.%s.bedgraph', mark_dir, version)
        write.table(coverage_mates[rCN %in% c('1|1', '1|2')][, c('chr', 'start', 'end', 'tCN')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)    
        return(coverage_mates)
    }

    iterative_gam = function(tab4reg) ## seems the second round will not add gain
    {
        k = 20; ki = 20
        bin_size = tab4reg$end[1] - tab4reg$start[1] + 1
        if(bin_size < 50E3) {k = 10; ki = 10}
        # my_gam <- gam(ceiling(bulk/tCN) ~ te(gcper, k=k) + te(mapscore, k=k) + te(cutSitesNumber, k=k) + te(comp_rank, k=k), data=tab4reg, family=nb(link=log), gamma=1E-5)    
        my_gam <- mgcv::gam(ceiling(bulk/tCN) ~ te(gcper, k=k) + te(cutSitesNumber, k=k) + te(comp_rank, k=k), data=tab4reg, family=nb(link=log), gamma=1E-5) ## mapscore seems to have too big variation, can remove it for regression

        tab4reg$resi_gam = trim_outlier(residuals(my_gam, type="working"))

        tab4reg$bulk_gam = ceiling(exp(predict(my_gam, newdata=tab4reg))*tab4reg$tCN)
        tab4reg$cnv_gam = tab4reg$bulk / tab4reg$bulk_gam * tab4reg$tCN
        tab4reg[, mean_gam_chunk:=my_mean(cnv_gam), by='chunk']
        res = list(tab4reg=tab4reg, my_gam=my_gam)
        return(res)
    }

    ## split_key: based on which column to split (eg: rCN)
    define_chunk = function(coverage_mates, maxgap=100E6, split_key) ## split each chr into chunks with the same rCN. maxgap: assign to different chunk if two same cnv_cat are too far away
    {
        chrs2look = unique(coverage_mates$chr)
        coverage_mates = foreach(chr2look=chrs2look, .combine=rbind) %do%
        {
            coverage_mates_chr = coverage_mates[chr==chr2look][order(start)]
            coverage_mates_chr$index = 1:nrow(coverage_mates_chr)
            chunk_label = paste0('g_', rleid(coverage_mates_chr[[split_key]]))
            coverage_mates_chr$chunk = chunk_label

            ################################

            ## connect CRC with same rCN if they are not far away from each other
            coverage_mates_chr_li = split(coverage_mates_chr, by=split_key)
            new_chunk = foreach(j=1:length(coverage_mates_chr_li), .combine=rbind) %do% 
            {
                tmp_by_chunk = split(coverage_mates_chr_li[[j]], by='chunk')
                start = sapply(tmp_by_chunk, function(v) min(v$index))
                end = sapply(tmp_by_chunk, function(v) max(v$index))
                GR = IRanges(start=start, end=end)
                GR_merged = reduce(GR, min.gapwidth=maxgap/bin_size2look)
                assign2chunk = subjectHits(findOverlaps(GR, GR_merged))
                for(z in 1:length(assign2chunk)) tmp_by_chunk[[z]]$chunk = sprintf('%s_c%s_%s', chr2look, j, assign2chunk[z])
                new_chunk = do.call(rbind, tmp_by_chunk)
                new_chunk
            }
            new_chunk
        }
        coverage_mates[, index:=NULL]

        ################################

        colors <- col2hex(c("red", "orange", "yellow", "green", "grey", "violet", 'yellow', "cyan", "magenta", "pink"))
        colors = rep(colors, 100)

        bedgraph = cbind(coverage_mates[, mget(c(GR_3_cols, 'chunk', 'seg_ratio'))], dot='.', coverage_mates[, mget(GR_3_cols[2:3])])
        bedgraph$col = colors[match(bedgraph$chunk, unique(bedgraph$chunk))]

        bed_f = sprintf('%s/cnv_chunk.v1.bed', mark_dir)
        write.table(bedgraph, file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
        return(coverage_mates)
    }

    regress_fun = function(coverage_mates, which2use, version)
    {
        # tab4reg = coverage_mates[tCN!=0][segLen > 5][rCN %in% which2use]

        tab4reg = coverage_mates[rCN %in% which2use][segLen > 5]
        tab_n4reg = coverage_mates[!tab4reg, on=GR_3_cols] 

        cat('Ratio4reg:\n')
        print(nrow(tab4reg)/nrow(coverage_mates))      

        iterative_gam_res = iterative_gam(tab4reg)
        tab4reg = iterative_gam_res$tab4reg
        my_gam = iterative_gam_res$my_gam

        cat('abs(tab4reg$resi_gam):\n')
        print(summary(abs(tab4reg$resi_gam)))

        bed_f = sprintf('%s/residual.%s.bedgraph', mark_dir, version)
        write.table(tab4reg[, c('chr', 'start', 'end', 'resi_gam')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

        coverage_mates$cnv_gam_del = ceiling(exp(predict(my_gam, newdata=coverage_mates))*coverage_mates$tCN) ## cnv_corrected gam. This value is very big, so should not worry about "ceiling"
        coverage_mates$cnv_gam = coverage_mates$bulk / coverage_mates$cnv_gam_del * coverage_mates$tCN
        coverage_mates[chr=='chrX']$cnv_gam = coverage_mates[chr=='chrX']$cnv_gam*factor4chrX
        coverage_mates[, cnv_gam_del:=NULL]

        bed_f = sprintf('%s/cnv_gam.%s.bedgraph', mark_dir, version)
        if(version=='v1') write.table(coverage_mates[rCN %in% c('1|1', '1|2')][, c('chr', 'start', 'end', 'cnv_gam')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
        if(version=='v2') write.table(coverage_mates[, c('chr', 'start', 'end', 'cnv_gam')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

        if(version=='v1') cat('generated cnv_gam for 11 12 region\n')
        if(version=='v2') cat('generated cnv_gam for 11 12 and all other regions\n')

        cols2keep = setdiff(colnames(coverage_mates), c('comp_rank', 'cutSitesNumber', 'gcper', 'mapscore'))
        if(version=='v2') coverage_mates = coverage_mates[, mget(cols2keep)]

        return(coverage_mates)
    }

    update_tCN = function(coverage_mates2use) ## update rCN, tCN by the regressed gam
    {
        cat('define the final rCN/aCN for all phased regions\n')
        get_cnv_range = function(cnv_gam) ## range to consider for a cnv_gam. Looks good by checking: x = seq(1, 10, by=0.3); cbind(x, do.call(rbind, lapply(x, get_cnv_range)))
        {
            # cnv_min = round_to_nearest(cnv_gam - ifelse(cnv_gam > 4, 1, 0.5))
            # cnv_max = round_to_nearest(cnv_gam + ifelse(cnv_gam > 4, 1, 0.5))  
            # cnv_min = floor(cnv_gam)
            # cnv_max = ceiling(cnv_gam)

            min_val = cnv_gam / 1.2 ## 20% difference
            max_val = cnv_gam * 1.2

            cnv_min = ceiling(min_val)
            cnv_max = floor(max_val)
            c(cnv_min=cnv_min, cnv_max=cnv_max)
        }

        get_feasible_rCN = function(cnv_gam, seg_ratio)
        {
            cnv_candidates = get_cnv_range(cnv_gam)['cnv_min'] : get_cnv_range(cnv_gam)['cnv_max']
            cnv_candidates = setdiff(cnv_candidates, 0)
            ks = foreach(k=1:nrow(hap_config_tab_relaxed), .combine=c) %do% ## if cnv_candidates = n_total x some integer
            {
                candidate_in = 1*(cnv_candidates %% hap_config_tab_relaxed$n_total[k] == 0)
                if(any(candidate_in==1)) return(k)
                return(0)
            }

            ks = setdiff(ks, 0)
            if(length(ks)==0) return(list(rCN=gsub('.mix.mix', '.mix', paste0(tab2use_i$rCN[1], '.mix')), tCN=cnv_gam))
          
            feasible_tab = hap_config_tab_relaxed[ks]
            feasible_tab = feasible_tab[abs(log(ratio) - log(seg_ratio)) <= log(1.1)]
            if(nrow(feasible_tab)==0) return(list(rCN=gsub('.mix.mix', '.mix', paste0(tab2use_i$rCN[1], '.mix')), tCN=cnv_gam))
            
            feasible_tab$seg_ratio = seg_ratio
            feasible_tab$cnv_gam = cnv_gam
            feasible_tab$deltaRatio = abs(seg_ratio - feasible_tab$ratio)

            feasible_tab = feasible_tab[!(n_mated=='1|2' & seg_ratio > 0.71)] ## this is a threshold not well defined by feasible_tab[abs(log(ratio) - log(seg_ratio)) <= log(1.1)], 24/09/2023
            if(nrow(feasible_tab)==0) return(list(rCN=gsub('.mix.mix', '.mix', paste0(tab2use_i$rCN[1], '.mix')), tCN=cnv_gam))
           
            feasible_tab$deltaCN = sapply(feasible_tab$n_total, function(v) min(abs(seq(1:20)*v - cnv_gam)))
            feasible_tab[, toChoose:=1*(deltaRatio==min(deltaRatio)), by='deltaCN'] ## sometimes, the same tCN can match to multiple aCNs, such as 4|5, 3|6
            feasible_tab = feasible_tab[toChoose==1]
            feasible_tab = feasible_tab[which.min(deltaCN)]

            tCN_tmp = feasible_tab$n_total*seq(1:20)
            tCN = tCN_tmp[which.min(abs(cnv_gam - tCN_tmp))]
            rCN = sprintf('%s|%s', feasible_tab$n_mat, feasible_tab$n_pat)
            out = list(rCN=rCN, tCN=tCN)
            return(out)
        }

        cols2keep = colnames(coverage_mates2use)
        tab_phased = coverage_mates2use[rCN!='n/p']
        tab_np = coverage_mates2use[rCN=='n/p']
        tab_np$tCN = tab_np$mean_gam_chunk
        li = split(tab_phased, by='chunk')

        tab_phased = foreach(i=1:length(li), .combine=rbind) %do%
        {
            tab2use_i = li[[i]]
            cnv_gam = tab2use_i$mean_gam_chunk[1]*1.1 ## some of the chromosomes, the coverage seems to be lower than expected, for example, chrX. So add this compensation / 24/09/2023
            seg_ratio = get_seg_ratio(tab2use_i)
            tab2use_i$seg_ratio = seg_ratio

            rtCN = get_feasible_rCN(cnv_gam, seg_ratio)
            tab2use_i$rCN = rtCN[['rCN']]
            tab2use_i$tCN = rtCN[['tCN']]
            return( tab2use_i )
        }

        coverage_mates2use = rbind(tab_phased[, mget(cols2keep)], tab_np[, mget(cols2keep)])
        
        bed_f = sprintf('%s/tCN.v3.bedgraph', mark_dir)
        write.table(coverage_mates2use[, c('chr', 'start', 'end', 'tCN')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)    

        return(coverage_mates2use)
    }

    ################################

    get_seg_ratio = function(coverage_tab)
    {
        coverage_tab = coverage_tab[All >= quantile(coverage_tab$All, 0.25)]
        max(median(coverage_tab$mat_inAll), median(coverage_tab$pat_inAll))
    }
    

    build_tab = function(wk_dir, region2exclude, bin_size2look)
    {
        coverage_mates = foreach(mate=c(mates, 'bulk')[1:3], .combine=merge_all) %do%
        {
            bed_f = sprintf('%s/coverage.%s.bedgraph', mark_dir, mate)
            tab = fread(bed_f)[, 1:4]
            colnames(tab) = c(GR_3_cols, mate)
            tab[chr %in% chr_chars2look[1:23]]
        }

        coverage_mates$All = (coverage_mates$mat + coverage_mates$pat)

        coverage_mates$mat_inAll = coverage_mates$mat / coverage_mates$All
        coverage_mates$pat_inAll = coverage_mates$pat / coverage_mates$All

        bulk_thresh = mean(coverage_mates[bulk > 0]$bulk) / 10
        coverage_mates = coverage_mates[bulk > bulk_thresh]

        coverage_mates = coverage_mates[!region2exclude, on=GR_3_cols]

        ################################ remove empty region and problem region

        coverage_mates = coverage_mates[!region2exclude, on=GR_3_cols]

        ################################ CRC, segmented

        CRC_f = sprintf('%s/diffCNV.bedgraph', mark_dir)
        CRC = fread(CRC_f)
        colnames(CRC) = c(GR_3_cols, 'CRC')

        # already excluded when computing CRC_f
        # CRC = CRC[!low_snp_region, on = GR_3_cols] ## some small regions can have false CRC segmentation, they are indeed low SNP region; this should be moved in the part to define diffCNV

        coverage_mates = merge(coverage_mates, CRC, all.x=TRUE)
        coverage_mates[is.na(CRC)]$CRC = -1 ## region without proper CRC
        coverage_mates = coverage_mates[order(chr, start)]
        cols2keep = colnames(coverage_mates)

        ################################

        coverage_mates_li = split(coverage_mates, by=c('chr', 'CRC')) ## split into region of same CRC value
        coverage_mates = foreach(i=1:length(coverage_mates_li), .combine=rbind) %do%
        {
            coverage_mates_i = coverage_mates_li[[i]]
            if(is.na(median(coverage_mates_i$pat_inAll, na.omit=TRUE)) | coverage_mates_i$CRC[1]==-1) ## if no phased contact info is available
            {
                coverage_mates_i$seg_ratio = -1
                coverage_mates_i$rCN = 'n/p'
                return(coverage_mates_i)
            }

            seg_ratio = get_seg_ratio(coverage_mates_i)
            coverage_mates_i$seg_ratio = seg_ratio
            coverage_mates_i$rCN = 'mix'

            close2n = hap_config_tab[between(coverage_mates_i$seg_ratio[1], ratio_min, ratio_max)]$index ## in which hap config category
            if(length(close2n) > 0) close2n = close2n[1] ## in case there is an overlap between hap_config_tab profile, select the one with smaller total cnv (hap_config_tab needs to be ordered by n_total)
            if(length(close2n)==0)
            {
                coverage_mates_i$rCN = 'mix'
                return(coverage_mates_i)    
            }

            coverage_mates_i$rCN = hap_config_tab$n_mated[close2n]
            coverage_mates_i
        }

        ################################ only assign meaningful value to 1|1 and 1|2; all other places set as psudo value 1

        coverage_mates$tCN = 1
        coverage_mates[rCN=='1|1']$tCN = 2
        coverage_mates[rCN=='1|2']$tCN = 3
        # coverage_mates[rCN=='1|3']$tCN = 4
        # coverage_mates[rCN=='1|4']$tCN = 5
        # coverage_mates[rCN=='2|3']$tCN = 5

        coverage_mates[rCN=='1|1']$seg_ratio = 1/2
        coverage_mates[rCN=='1|2']$seg_ratio = 2/3

        coverage_mates[, segLen:=length(start)*bin_size2look/1E6, by=c('chr', 'seg_ratio')] ## combine segments in the same chr together

        ################################

        bedgraph = cbind(coverage_mates[, mget(c(GR_3_cols, 'rCN', 'seg_ratio'))], dot='.', coverage_mates[, mget(GR_3_cols[2:3])])
        bedgraph$col = hap_config_tab$col[match(bedgraph$rCN, hap_config_tab$n_mated)]

        for(i in 1:ncol(bedgraph)) bedgraph[[i]] = as.character(bedgraph[[i]])
        bed_f = sprintf('%s/CRC_cat.init.bed', mark_dir)
        fwrite(bedgraph, file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

        ################################

        coverage_mates[, CRC:=NULL]
        coverage_mates[, mat:=NULL]
        coverage_mates[, pat:=NULL]
       
        return(coverage_mates)
    }

    ################################ segmentation inside each chunk

    redefine_chunk = function(coverage_mates, which2use)
    {
        if(which2use[1]=='All') which2use = unique(coverage_mates$rCN)
        cat('n_chunk before:', length(unique(coverage_mates[rCN %in% which2use]$chunk)), '\n')
        interative_merge = function(coverage, seg_type) ## merge iteratively of neighboring segs
        {
            uniq_val = log2(unique(coverage$mean_gam_chunk))
            n = length(uniq_val) - 1
            if(n==0) return(coverage)
            abs_dev = abs(diff(uniq_val))
            min_cnv = sapply(1:n, function(v) min(uniq_val[v:(v+1)]))
            dev_thresh = log2(1.5) ## at least 1.5x difference in the same diffCNV region
            if(seg_type=='n/p') dev_thresh = log2(1.25) ## for unphased region, lower threshold. This might not be the best threshold to differ 4 - 5, 3 - 4

            dev_tab = data.table(min_cnv=min_cnv, abs_dev=abs_dev, dev_thresh=dev_thresh, index=1:n)
            dev_tab = dev_tab[abs_dev < dev_thresh] ## to merge

            if(nrow(dev_tab)==0) return(coverage)

            index2merge = dev_tab$index[which.min(dev_tab$abs_dev)]
            index2merge = c(index2merge, index2merge + 1)
            coverage2merge = coverage[log2(mean_gam_chunk) %in% uniq_val[index2merge]]

            coverage_n2merge = coverage[!coverage2merge, on=GR_3_cols]
            coverage2merge[, mean_gam_chunk:=my_mean(cnv_gam)]
            coverage = rbind(coverage2merge, coverage_n2merge)[order(start)]
            return(interative_merge(coverage, seg_type))
        } 
        
        coverage_mates_li = split(coverage_mates, by='chunk')
        coverage_mates = foreach(i=1:length(coverage_mates_li), .combine=rbind) %do%
        {
            seg_type = coverage_mates_li[[i]]$rCN[1]
            coverage = coverage_mates_li[[i]][order(start)]
            coverage$index = 1:nrow(coverage)
            cnv_tree <- rpart(cnv_gam ~ index, data=coverage, control=rpart.control(minsplit=2, minbucket=2, cp=0.01)) ## segmentation
            coverage$seg_val = predict(cnv_tree, data.table(index=coverage$index)) ## different seg_val means different segmentation
            coverage[, mean_gam_chunk:=my_mean(cnv_gam), by='seg_val'] ## instead of using the original segmentation value, update here to avoid some super high values dominate the seg value
            coverage = interative_merge(coverage, seg_type)
            coverage$chunk = paste0(coverage$chunk, '.', as.integer(factor(coverage$mean_gam_chunk)))
            coverage[, index:=NULL]
            coverage
        }

        colors <- col2hex(c("red", "orange", "yellow", "green", "grey", "violet", 'yellow', "cyan", "magenta", "pink"))
        colors = rep(colors, 100)
        bedgraph = cbind(coverage_mates[, mget(c(GR_3_cols, 'chunk', 'seg_ratio'))], dot='.', coverage_mates[, mget(GR_3_cols[2:3])])
        bedgraph$col = colors[match(bedgraph$chunk, unique(bedgraph$chunk))]

        bed_f = sprintf('%s/cnv_chunk.v2.bed', mark_dir)
        write.table(bedgraph, file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)        

        cat('n_chunk after:', length(unique(coverage_mates[rCN %in% which2use]$chunk)), '\n')
        return(coverage_mates)     
    }

    ################################

    segmentation_fun = function(coverage_mates)
    {
        interative_merge = function(coverage, seg_type) ## merge iteratively of neighboring segs
        {
            dev_thresh = ifelse(seg_type %in% c('n/p', 'mix'), 0.1, 0.15) ## deviation thresh from the mean of the two segments; higher thresh for phased homogenious seg

            uniq_val = unique(coverage$mean_gam_chunk)
            n = length(uniq_val) - 1
            # cat('n_seg:', n + 1, '\n')
            if(n==0) return(coverage)
            deviations = sapply(1:n, function(v)
            {
                x = uniq_val[v]
                y = uniq_val[v+1]
                mean_xy = (x + y) /2
                out = (max(x, y) - mean_xy) / mean_xy
                return(out)
            }) ## deviation from the center

            min_dev = min(deviations)
            if(min_dev > dev_thresh) return(coverage)

            index2merge = which(deviations==min_dev)
            index2merge = c(index2merge, index2merge + 1)
            coverage2merge = coverage[mean_gam_chunk %in% uniq_val[index2merge]]

            coverage_n2merge = coverage[!coverage2merge, on=GR_3_cols]
            coverage2merge[, mean_gam_chunk:=my_mean(cnv_gam)]
            coverage = rbind(coverage2merge, coverage_n2merge)[order(start)]
            return(interative_merge(coverage, seg_type))
        } 
        

        coverage_mates_li = split(coverage_mates, by='chunk')
        coverage_mates = foreach(i=1:length(coverage_mates_li), .combine=rbind) %do%
        {
            seg_type = coverage_mates_li[[i]]$rCN[1]
            coverage = coverage_mates_li[[i]][order(start)]
            coverage$index = 1:nrow(coverage)
            cnv_tree <- rpart(cnv_gam ~ index, data=coverage, control=rpart.control(minsplit=5, minbucket=5, cp=0.01)) ## segmentation
            coverage$seg_val = predict(cnv_tree, data.table(index=coverage$index)) ## different seg_val means different segmentation
            coverage[, mean_gam_chunk:=my_mean(cnv_gam), by='seg_val'] ## instead of using the original segmentation value, update here to avoid some super high values dominate the seg value
            coverage = interative_merge(coverage, seg_type)
            coverage$chunk = paste0(coverage$chunk, '.', as.integer(factor(coverage$mean_gam_chunk)))
            coverage[, index:=NULL]
            coverage
        }
        return(coverage_mates)     
    }

    ################################

    define_aCN = function(coverage_mates2use)
    {
        cat('generate aCN.mat and aCN.pat for phased non-mix region\n')
        coverage_mates2use$aCN = 'TBD' ## cnv label
        coverage_mates2use_li = split(coverage_mates2use, by=c('rCN', 'tCN', 'chunk'))

        for(i in 1:length(coverage_mates2use_li))
        {
            coverage_mates2use_i = coverage_mates2use_li[[i]]
            rCN = coverage_mates2use_i$rCN[1]
            if(rCN=='n/p') 
            {
                coverage_mates2use_li[[i]]$aCN = 'n/p'
                next
            }

            tCN = coverage_mates2use_i$tCN[1]
            if(rCN=='mix.mix') cnv_mates_num = numeric()
            if(rCN!='mix.mix') cnv_mates_num = sort(as.numeric(strsplit(gsub('.mix', '|', rCN), '[|]')[[1]]))

            is_pat_more = median(coverage_mates2use_i$mat_inAll) <= median(coverage_mates2use_i$pat_inAll)
            if(!is_pat_more) cnv_mates_num = rev(cnv_mates_num)

            cnv_mates_num = tCN / sum(cnv_mates_num) * cnv_mates_num
            rCN_new = paste0(cnv_mates_num, collapse='|')
            if(grepl('mix', rCN)) rCN_new = paste0(paste0(cnv_mates_num, collapse=''), '.mix')
            coverage_mates2use_li[[i]]$aCN = rCN_new
        }

        coverage_mates2use = do.call(rbind, coverage_mates2use_li)
        coverage_mates2use[grepl('mix', rCN)]$aCN = 'mix'
        coverage_mates2use[grepl('mix', aCN)]$aCN = 'mix'

        return(coverage_mates2use)
    }

    ################################

    define4LOH = function(coverage_mates_final)
    {
        coverage_mates_final$index = 1:nrow(coverage_mates_final)
        LOH_region = coverage_mates_final[LOH==1]
        # LOH_region[, seg_ratio:=mean(pat_inAll, na.omit=TRUE), by='LOH']

        LOH_region[, direction_by_chr:=median(mat_inAll, na.rm=TRUE) < median(pat_inAll, na.rm=TRUE), by='chr'] ## assign a 'direction' at chromosome level, so that to keep the same direction for unphased LOH
        dir_by_chr = unique(LOH_region[, c('chr', 'direction_by_chr')])
        dir_by_chr[is.na(dir_by_chr)] = 1

        LOH_region_li = split(LOH_region, by=c('chunk', 'LOH'))
        LOH_region = foreach(i=1:length(LOH_region_li), .combine=rbind) %do%
        {
            LOH_region_i = LOH_region_li[[i]]
            direction = dir_by_chr[chr==LOH_region_i$chr[1]]$direction_by_chr
            mean_cnv_gam = mean(LOH_region_i$cnv_gam)
            nearest = round_to_nearest(mean_cnv_gam)
            mean_cnv_gam = ifelse(my_dev_fun(nearest, mean_cnv_gam) < log(1.2), nearest, mean_cnv_gam)
            aCN = substr(as.character(mean_cnv_gam), 1, 3)
            # direction = median(LOH_region_i$mat_inAll, na.rm=TRUE) < median(LOH_region_i$pat_inAll, na.rm=TRUE)
            # if(is.na(direction)) direction = 1 ## assign to pat if no any phased reads info
            if(direction) aCN = sprintf('0|%s', aCN)
            if(!direction) aCN = sprintf('%s|0', aCN)

            LOH_region_i$aCN = aCN
            LOH_region_i
        }

        LOH_region[, direction_by_chr:=NULL]

        coverage_mates_final = rbind(coverage_mates_final[index %in% setdiff(coverage_mates_final$index, LOH_region$index)], LOH_region)
        coverage_mates_final[grepl('[|]', aCN)]$tCN = apply(do.call(rbind, strsplit(coverage_mates_final[grepl('[|]', aCN)]$aCN, '[|]')), 1, function(v) sum(as.numeric(v)))

        coverage_mates_phased = coverage_mates_final[grepl('[|]', aCN)]
        coverage_mates_phased$cnv_mat = as.numeric(do.call(rbind, strsplit(coverage_mates_phased$aCN, split='[|]'))[,1])
        coverage_mates_phased$cnv_pat = as.numeric(do.call(rbind, strsplit(coverage_mates_phased$aCN, split='[|]'))[,2])

        coverage_mates_phased = foreach(chr2look=unique(coverage_mates_phased$chr), .combine=rbind) %do%
        {
            coverage_mates_phased_chr = coverage_mates_phased[chr==chr2look]
            coverage_mates_phased_chr$cnv_mat.tmp = coverage_mates_phased_chr$cnv_mat
            coverage_mates_phased_chr$cnv_pat.tmp = coverage_mates_phased_chr$cnv_pat
            if(sum(coverage_mates_phased_chr$cnv_mat) > sum(coverage_mates_phased_chr$cnv_pat)) setnames(coverage_mates_phased_chr, c('cnv_mat.tmp', 'cnv_pat.tmp'), c('cnv_pat.tmp', 'cnv_mat.tmp')) 
            coverage_mates_phased_chr
        }


        bed_f = sprintf('%s/aCN.mat.bedgraph', mark_dir)       
        write.table(coverage_mates_phased[, c('chr', 'start', 'end', 'cnv_mat')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)    
        bed_f = sprintf('%s/aCN.pat.bedgraph', mark_dir)
        write.table(coverage_mates_phased[, c('chr', 'start', 'end', 'cnv_pat')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)    

        bed_f = sprintf('%s/aCN.mat4igv.bedgraph', mark_dir)       
        write.table(coverage_mates_phased[, c('chr', 'start', 'end', 'cnv_mat.tmp')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)    
        bed_f = sprintf('%s/aCN.pat4igv.bedgraph', mark_dir)
        write.table(coverage_mates_phased[, c('chr', 'start', 'end', 'cnv_pat.tmp')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)    
   
        return(coverage_mates_final)
    }


 