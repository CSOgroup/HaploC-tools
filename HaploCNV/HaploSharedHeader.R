   compute_coverage_sparsity = function(wk_dir, map_type, bin_size=100E3)
    {
        ################################ genarete coverage and sparsity

        if(map_type=="phased") which_map = mates
        if(map_type!="phased") which_map = 'bulk'

        for(mate in which_map) get_hic_coverage(wk_dir, bin_size=bin_size, mate, flanking_bp=NA, impute=TRUE)   
        for(mate in which_map) get_hic_sparsity(wk_dir, bin_size=bin_size, mate)
        
        ################################

        if(map_type=='phased')
        {
            coverage_tab_mates = foreach(mate=mates) %do% 
            {
                coverage_f = file.path(wk_dir, sprintf("CNV/%skb/coverage.%s.bedgraph", bin_size/1E3, mate))
                coverage_tab = fread(coverage_f)[, 1:4]
                colnames(coverage_tab) = c("chr", "pos_start", "pos_end", paste0("coverage_", mate))
                coverage_tab
            }

            coverage_tab_mate = merge(coverage_tab_mates[[1]], coverage_tab_mates[[2]], by=c("chr", "pos_start", "pos_end"), all=TRUE)
            coverage_tab_mate$total = coverage_tab_mate$coverage_mat + coverage_tab_mate$coverage_pat
            coverage_tab_mate = coverage_tab_mate[total > 0]
            coverage_tab_mate$ratio_mat = coverage_tab_mate$coverage_mat / coverage_tab_mate$total - 0.5
            coverage_tab_mate$ratio_pat = coverage_tab_mate$coverage_pat / coverage_tab_mate$total - 0.5
            
            no_out = foreach(mate=which_map) %do% 
            {
                ratio_f = file.path(wk_dir, sprintf("CNV/%skb/coverage.ratio.%s.bedgraph", bin_size/1E3, mate))
                cols2keep = c("chr", "pos_start", "pos_end", paste0("ratio_", mate))
                fwrite(coverage_tab_mate[, mget(cols2keep)], file=ratio_f, sep='\t', quote=FALSE, col.names=FALSE)
            }
        }
    }

    get_hic_sparsity = function(wk_dir, bin_size=100E3, mate) ## get contact mat sparsity
    {
        if(mate!='bulk') hic_f = get_mate_hic(wk_dir, mate)
        if(mate=='bulk') hic_f = get_bulk_hic(wk_dir)
        chr_nums2use = get_chr2use(wk_dir)

        coverage_chrs = foreach(chr_num2look=chr_nums2use, .combine=rbind) %do%
        {
            if(chr_num2look==23) chr_num2look = 'X'
            chr_char = paste0('chr', chr_num2look)
            coverage_bed = foreach(mate=mate, .combine=rbind) %do%
            {   
                chr2strawr = as.character(chr_num2look)
                dump_mat = try(silent=TRUE, strawr::straw(matrix='observed', binsize=bin_size, chr1loc=chr2strawr, chr2loc=chr2strawr, unit='BP', fname=hic_f, norm='NONE'))
                if(class(dump_mat)[1]=="try-error") return(NULL)
                
                dump_mat = as.data.table(dump_mat)
                dump_mat = na.omit(dump_mat)
                colnames(dump_mat) = c("pos_1", "pos_2", "counts")
                dump_mat$counts = 1 ## for sparsity 

                mat = convert_contact_format_fun(dump_mat, from='format_dump', to='format_sparse', chr_num=chr_num2look, binsize=bin_size)$format_sparse

                coverage = Matrix::rowSums(mat) / nrow(mat) ## sparsity
                bin_index = as.numeric(names(coverage))
                coverage_tab = data.table(chr=chr_char, pos_start=as.character((bin_index-1)*bin_size + 1), pos_end=as.character((bin_index)*bin_size), val=unname(coverage), mate=mate)
                return(coverage_tab)
            }
            coverage_bed
        }
        sparsity_f = file.path(wk_dir, sprintf("CNV/%skb/sparsity.%s.bedgraph", bin_size/1E3, mate))
        fwrite(coverage_chrs, file=sparsity_f, sep='\t', quote=FALSE, col.names=FALSE)
    }


    get_hic_coverage = function(wk_dir, bin_size, mate, flanking_bp=NA, type='hic', impute=TRUE)
    {
        if(mate=='bulk') hic_f = get_bulk_hic(wk_dir)
        if(mate!='bulk') hic_f = get_mate_hic(wk_dir, mate, impute=impute)

        chr_nums2use = get_chr2use(wk_dir)
        coverage_chrs = foreach(chr_num2look=chr_nums2use, .combine=rbind) %do%
        {
            if(chr_num2look==23) chr_num2look = 'X'
            chr_char = paste0('chr', chr_num2look)
            coverage_bed = foreach(mate=mate, .combine=rbind) %do%
            {
                chr2strawr = as.character(chr_num2look)
                dump_mat = try(silent=TRUE, strawr::straw(matrix='observed', binsize=bin_size, chr1loc=chr2strawr, chr2loc=chr2strawr, unit='BP', fname=hic_f, norm='NONE'))
                if(class(dump_mat)[1]=="try-error") return(NULL)
                dump_mat = as.data.table(dump_mat)
                dump_mat = na.omit(dump_mat[abs(x - y) <= 10E6])
                colnames(dump_mat) = c("pos_1", "pos_2", "counts")
                
                mat = convert_contact_format_fun(dump_mat, from='format_dump', to='format_sparse', chr_num=chr_num2look, binsize=bin_size)$format_sparse
                # mat[abs(row(mat) - col(mat)) > (10E6/bin_size)] = NA ## only keep local contacts / before 04.05.2023 / keep all within the same flank region 04.05.2023, this can reduce noise 
                # coverage = apply(mat, 1, sum, na.rm=TRUE)
               
                coverage = Matrix::rowSums(mat)
                coverage[coverage > (quantile(coverage, 0.8)*100)] = (quantile(coverage, 0.8)*100) ## filter the extreme values

                bin_index = as.numeric(names(coverage))
                coverage_tab = data.table(chr=chr_char, pos_start=as.character((bin_index-1)*bin_size + 1), pos_end=as.character((bin_index)*bin_size), val=unname(coverage), mate=mate)
                return(coverage_tab)
            }
            coverage_bed
        }
  
        coverage_f = file.path(wk_dir, sprintf("CNV/%skb/coverage.%s.bedgraph", bin_size/1E3, mate))
        fwrite(coverage_chrs, file=coverage_f, sep='\t', quote=FALSE, col.names=FALSE)
        return(coverage_chrs)
    }



