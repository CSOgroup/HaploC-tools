    options(stringsAsFactors=FALSE)
       
    ##################################### packages
       
    try(silent=TRUE, suppressMessages(library(data.table)))
    try(silent=TRUE, suppressMessages(library(doParallel)))
    try(silent=TRUE, suppressMessages(library(GenomicRanges)))
    try(silent=TRUE, suppressMessages(library(ggplot2)))
    # try(silent=TRUE, suppressMessages(library(ggpubr)))
    try(silent=TRUE, suppressMessages(library(gplots))) 
    try(silent=TRUE, suppressMessages(library(grid)))
    try(silent=TRUE, suppressMessages(library(Matrix)))
    try(silent=TRUE, suppressMessages(library(mgcv)))
    try(silent=TRUE, suppressMessages(library(rpart)))
    try(silent=TRUE, suppressMessages(library(strawr)))
    try(silent=TRUE, suppressMessages(library(tools)))
    try(silent=TRUE, suppressMessages(library(viridis)))

    ##################################### dir

    HaploC_dir = file.path(repo_dir, 'HaploC')
    HaploCNV_dir = file.path(repo_dir, 'HaploCNV')
    juicerDir = file.path(repo_dir, 'HaploC/dependancies/juicer')
    diffComp_dir = file.path(repo_dir, 'diffComp')
    diffIns_dir = file.path(repo_dir, 'diffIns')
    downstream_dir= file.path(repo_dir, 'downstream_analysis')

    ##################################### consts

    ref_fa.hg19 = sprintf("%s/genomicData/references/hg19_chr1_22XYM.fa", repo_dir)
    ref_fa.mm10 = sprintf("%s/genomicData/references/mm10.fa", repo_dir)
    ref_fa_li = list(hg19=ref_fa.hg19, mm10=ref_fa.mm10)

    chrom_sizes = file.path(repo_dir, 'genomicData/references/hg19.chrom.sizes')
    mates = c("mat", "pat"); GR_3_cols = chr_start_end = c('chr', 'start', 'end'); chr_bin = c('chr', 'bin_index')
    comp_names = c("A.1.1", "A.1.2", "A.2.1", "A.2.2", "B.1.1", "B.1.2", "B.2.1", "B.2.2"); comp_vals = (8:1)/8; names(comp_vals) = comp_names
    comp_cols = c("#FF0000", "#FF4848", "#FF9191", "#FFDADA", "#DADAFF", "#9191FF", "#4848FF", "#0000FF"); names(comp_cols) = comp_names

    ##################################### excutable

    Rscript = 'Rscript'

    ##################################### functions

    my_hclust = function(...) rev(as.dendrogram(hclust(..., method='average')))

    my_forceSymmetric = function(mat) ## Matrix package does not work on UNIL HPC conda env intg
    {
        new_mat = mat + t(mat)
        diag(new_mat) = diag(mat)
        return(new_mat)
    }

    get_bulk_hic = function(wk_dir)
    {
        fs = list.files(sprintf('%s/mega/aligned', wk_dir), full.names=TRUE)
        hic_f = fs[substr(fs, nchar(fs)-2, nchar(fs))=='hic']
        return(hic_f)
    }

    get_mate_hic = function(wk_dir, mate, impute=TRUE)
    {
        hic_f = sprintf('%s/mega/HapCut/hic_phased/all_chrs.opt.%s.hic', wk_dir, mate)
        if(!impute) hic_f = sprintf('%s/mega/HapCut/hic_phased/all_chrs.opt.n_impute.%s.hic', wk_dir, mate)
        return(hic_f)
    }

    generate_region_tab = function(bin_size)
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

    straw_from_hic = function(hic_f, chr_num, bin_size, norm, ob_oe) ## '1', ..., 'X'
    {
        chr_name = paste0("chr", chr_num)
        chr2strawr = chr_num
        if(any(grepl('chr', readHicChroms(hic_f)$name))) chr2strawr = chr_name

        if(norm=='norm')
        {
            mat_dump = try(silent=TRUE, strawr::straw(matrix=ob_oe, binsize=bin_size, chr1loc=chr_num, chr2loc=chr_num, unit='BP', fname=hic_f, norm='KR') %>% data.table)
            if(class(mat_dump)=='try-error') mat_dump = try(silent=TRUE, strawr::straw(matrix=ob_oe, binsize=bin_size, chr1loc=chr_num, chr2loc=chr_num, unit='BP', fname=hic_f, norm='VC_SQRT') %>% data.table)
            if(class(mat_dump)=='try-error') return(NULL)
            colnames(mat_dump) = c("pos_1", "pos_2", "counts")
            mat_dump = mat_dump[!is.nan(counts)]
        }

        if(norm!='norm')
        {
            mat_dump = try(silent=TRUE, strawr::straw(matrix=ob_oe, binsize=bin_size, chr1loc=chr_num, chr2loc=chr_num, unit='BP', fname=hic_f, norm='NONE') %>% data.table)
            if(class(mat_dump)=='try-error') return(NULL)
            colnames(mat_dump) = c("pos_1", "pos_2", "counts")
            mat_dump = mat_dump[!is.nan(counts)]
        }

        if(nrow(na.omit(mat_dump)) < 100) mat_dump = NULL
        return(mat_dump)
    }

    get_chr2use = function(wk_dir, genome_version=NULL) ## do not analyze chrs with less than 1E3 SNPs // tmp
    {
        chr2phase_f = file.path(wk_dir, 'mega/VCFs/chr2phase.txt')
        if(!file.exists(chr2phase_f))
        {
            snp_count_f = file.path(wk_dir, 'mega/VCFs/n_variants_het.no_ins.txt')
            snp_count_tab = read.table(snp_count_f)
            colnames(snp_count_tab) = c('chr', 'n', 'QUAL')
            chr_num2use = as.character(unique(subset(snp_count_tab, QUAL == 10 & n > 1E3)$chr)) ## IMR90 chr22 is 5E4
            if(genome_version=='hg19') chr_num2use[chr_num2use=='23'] = 'X'
            if(genome_version=='mm10') chr_num2use[chr_num2use=='20'] = 'X'            
            cat(chr_num2use, file=chr2phase_f)
        }
        chr_num2use = unname(strsplit(readLines(chr2phase_f, warn=FALSE), ' ')[[1]])
        return(chr_num2use)
    }


    convert_contact_format_fun = function(input_dat=NULL, chr_num, from='format_dump', to='format_sparse', binsize=20E3)
    {
        chr_name = paste0('chr', chr_num)
        contact_formats = list()

        if(from=='format_dump') 
        {
         if( !all(colnames(input_dat)==c('pos_1', 'pos_2', 'counts')) ) stop('Your input_dat should have colnames as: pos_1, pos_2, counts')
         contact_formats[['format_dump']] = input_dat
        }

        if( from=="format_dump" & to=="format_pre" )
        {
         format_pre = data.table(str1=0, chr1=chr_num, pos_1=as.character(input_dat$pos_1), frag1=0, str2=0, chr2=chr_num, pos_2=as.character(input_dat$pos_2), frag2=1, score=input_dat$counts)
         format_pre = format_pre[order(format_pre$pos_1, format_pre$pos_2), ]
         contact_formats[['format_pre']] = format_pre
         return(contact_formats)
        }

        #########################################

        ## my old script to convert dump format to sparse_mat format
        if(to=='format_sparse')
        {
         combined_xk_oe_raw = contact_formats$format_dump
         combined_xk_oe_raw[,1] = combined_xk_oe_raw[,1]/binsize
         combined_xk_oe_raw[,2] = combined_xk_oe_raw[,2]/binsize
         combined_xk_oe = combined_xk_oe_raw

         colnames(combined_xk_oe) = c('pos_1', 'pos_2', 'val') 
         if(!all(combined_xk_oe[[2]] >= combined_xk_oe[[1]])) stop('\nYou provided matrix does not represent an upper triangular matrix!\n\n')

         oe_size = max(max(combined_xk_oe[[1]]), max(combined_xk_oe[[2]])) + 1 ## should +1 because combined_xk_oe index starts from 0 (bin 0 represents: 0-10E3, checked by looking at the juicebox map, 2018-11-19)
         mat_oe_sparse = Matrix(0, nrow=oe_size, ncol=oe_size)
         mat_oe_sparse[cbind(combined_xk_oe[[1]]+1, combined_xk_oe[[2]]+1)] <- combined_xk_oe[[3]]

         rownames(mat_oe_sparse) = colnames(mat_oe_sparse) = as.character( 1:nrow(mat_oe_sparse) )
         # mat_oe_sparse = my_forceSymmetric(mat_oe_sparse, uplo='U')
         mat_oe_sparse = my_forceSymmetric(mat_oe_sparse)         
         contact_formats[['format_sparse']] = mat_oe_sparse
        }


        if(from=='format_sparse' & to=='format_dump')
        {
         contact_formats[['format_sparse']] = input_dat
         format_three_column = as.data.table(summary(contact_formats[['format_sparse']]))
         format_dump = data.table(binI=(format_three_column[[1]]-1)*binsize, binJ=(format_three_column[[2]]-1)*binsize, counts=format_three_column[[3]])
        }

        ######################################### dump to pre
        
        if( to=='format_pre' )
        {
         pos_I = do.call(rbind, strsplit(contact_formats$format_hicdc[['binI']], '-'))
         pos_J = do.call(rbind, strsplit(contact_formats$format_hicdc[['binJ']], '-'))
         values = contact_formats[['format_hicdc']]$counts

         format_pre = data.table(str1=0, chr1=chr_num, pos_1=as.character(as.numeric(pos_I[,3]) - binsize/2), frag1=0, str2=0, chr2=chr_num, pos_2=as.character(as.numeric(pos_J[,3]) - binsize/2), frag2=1, score=values)
         format_pre = format_pre[order(format_pre$pos_1, format_pre$pos_2), ]
         contact_formats[['format_pre']] = format_pre
        }

        return(contact_formats)
    }

    merge_all = function(...) merge(..., all=TRUE)
    