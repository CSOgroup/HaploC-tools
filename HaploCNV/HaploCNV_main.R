    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE) 
    wk_dir = args[1]
    bin_size2look = as.numeric(args[2])

    print(bin_size2look)

    # wk_dir = '/scratch/yliu5/1.Projects/15.Phasing/1.Data/WSU_v2.chr4_8_14'

    output_dir = sprintf('%s/CNV/%skb', wk_dir, bin_size2look/1E3)
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    mark_dir = sprintf('%s/CNV/%skb', wk_dir, bin_size2look/1E3)
   
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
    cutSitesNumber_gc_mapscore_rds_f = file.path(HaploCNV_dir, 'cutSitesNumber_gc_mapscore.rds')

    ################################

    # source(sprintf('%s/header/header.R', repo_dir))
    source(sprintf('%s/header_v2.R', HaploCNV_dir))
    source(sprintf('%s/helper.R', HaploCNV_dir))
    # source(sprintf('%s/HaploSharedHeader.R', HaploCNV_dir))     
    # source(sprintf('%s/diffComp_fun.R', diffComp_dir))
    # source(sprintf('%s/CALDER_fun.R', diffComp_dir))

    ################################

    chr_nums2look = get_chr2use(wk_dir)
    chr_chars2look = sprintf('chr%s', chr_nums2look)

    factor4chrX = 1.10 ## chrX seems have less contacts. Multiply this factor to raise it up

    ################################

    # for(map_type in c("phased", 'unphased')) compute_coverage_sparsity(ID, map_type, bin_size=bin_size2look)
    low_snp_region = get_low_snp_region(wk_dir, redo_count=TRUE)
    low_snp_region = expand_comp_tab_v2(low_snp_region, bin_size=bin_size2look); low_snp_region[, bin_index:=NULL]
    problem_region = get_encode_black_bins(bin_size=bin_size2look) ## problematic region to exclude
    # saveRDS(problem_region, file=file.path(post_script_dir, '7.Post_analysis/5.CNV/HaploCNV/ENCODE_blacklist/blacklist_bins.rds'))
    del_region = get_deletions(wk_dir)[, mget(GR_3_cols)] ## deleted region / no bulk contacts
    del_region = generate_binned_tab(makeGRangesFromDataFrame(del_region), bin_size2look)

    ################################

    region2exclude = unique(rbind(problem_region, del_region))

    ################################

    generate_diffCNV_seg(wk_dir, bin_size=bin_size2look, region2exclude=region2exclude, low_snp_region=low_snp_region) ## use 100kb to define diffCNV to reduce noise

    ################################

    hap_config_tab = get_hap_config_tab()
    hap_config_tab_relaxed = get_hap_config_tab_relaxed()
 
    ################################ get cnv cat freq. Indeed, 11 and 12 are most frequently represented

    coverage_mates = build_tab(wk_dir, region2exclude, bin_size2look) ## create coverage_mates profile
    coverage_mates = add_comp(wk_dir)
    coverage_mates = add_tech_features(coverage_mates, bin_size=bin_size2look, extream=0) ## all tCN got values. tCN=1: if rCN in 'mixed', or 'n/p'
    coverage_mates = define_chunk(coverage_mates, split_key='rCN', maxgap=10E6) ## define chunks based on rCN, and merge nearby same chunks

    ################################

    coverage_mates = set_aCN_11_12(coverage_mates, which2use='bulk', version='v1') ## decide tCN in 11, 22

    ################################ regression

    coverage_mates = regress_fun(coverage_mates, which2use=c('1|1', '1|2'), version='v1')

    ################################ regression

    ## repeat the process using corrected coverage
    ## check whether there is still different cnv in the same chunk
    ## the goal is not to define precisely the chunk but for better set_aCN_11_12
    coverage_mates = redefine_chunk(coverage_mates, which2use=c('1|1', '1|2'))
   
    coverage_mates = set_aCN_11_12(coverage_mates, which2use='cnv_gam', version='v2') ## this seems will not have an effect after tCN_11_12
    coverage_mates = regress_fun(coverage_mates, which2use=c('1|1', '1|2'), version='v2') ## abs(tab4reg$resi_gam) will become lower  
    coverage_mates2use = redefine_chunk(coverage_mates, which2use='All') ## cnv_gam has been properly defined
 
    ################################ 
  
    ## define tCN.v3
    coverage_mates2use = update_tCN(coverage_mates2use) ## update again cnv assignment, for phased region
    coverage_mates_final = define_aCN(coverage_mates2use)
    coverage_mates_final = add_LOH_info(coverage_mates_final, wk_dir, bin_size=bin_size2look)
    coverage_mates_final = define4LOH(coverage_mates_final)    

    ################################

    coverage_mates = coverage_mates_final
    coverage_mates.rCNed = coverage_mates[grepl('[|]', aCN)] ## rCN infered
    coverage_mates.cnv_n_guessed = coverage_mates[!grepl('[|]', aCN)] ## rCN not infered


    if(nrow(coverage_mates.cnv_n_guessed) > 0) ## assign an integer cnv to not phased region, if possible
    {
        coverage_mates.cnv_n_guessed$close2n = round_to_nearest(coverage_mates.cnv_n_guessed$mean_gam_chunk)
        coverage_mates.cnv_n_guessed$deviation = apply(coverage_mates.cnv_n_guessed[, c('close2n', 'mean_gam_chunk')], 1, function(v) my_dev_fun(v[1], v[2]))
        coverage_mates.cnv_n_guessed[deviation < log(1.2)]$tCN = coverage_mates.cnv_n_guessed[deviation < log(1.2)]$close2n
        coverage_mates.cnv_n_guessed[deviation >= log(1.2)]$tCN = coverage_mates.cnv_n_guessed[deviation >= log(1.2)]$mean_gam_chunk
        coverage_mates.cnv_n_guessed = coverage_mates.cnv_n_guessed[, mget(colnames(coverage_mates))]
    }

    coverage_mates = rbind(coverage_mates.rCNed, coverage_mates.cnv_n_guessed)

    ################################

    bed_f = sprintf('%s/tCN.bedgraph', mark_dir)
    rds_f = sprintf('%s/HaploCNV.rds', mark_dir)    
    write.table(coverage_mates[, c('chr', 'start', 'end', 'tCN')], file=bed_f, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
    saveRDS(coverage_mates, file=rds_f)

    ################################

