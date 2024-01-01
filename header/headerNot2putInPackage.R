    deep_clean_wk_dir_and_move2nas = function(LINE, wk_dir=NA, clean=FALSE, mv2nas=FALSE)
    {
        if(is.na(wk_dir)) wk_dir = file.path(phasing_dir, LINE)
        system(sprintf('du -sh %s', wk_dir))
        if(clean==FALSE) return("")
        fs_delete = list.files(wk_dir, recursive=TRUE, full.names=TRUE, pattern='delete')
        unlink(fs_delete)

        system(sprintf("rm -fr %s/*fastq.bp.gz", wk_dir))
        system(sprintf("rm -fr %s/chunk* -r", wk_dir))
        system(sprintf("rm -fr %s/mega/splits", wk_dir))
        system(sprintf("rm -fr %s/mega/aligned/merged_nodups.txt*", wk_dir))

        system(sprintf("rm -fr %s/mega/HapCut/fragment*", wk_dir))
        system(sprintf("rm -fr %s/mega/HapCut/shapeit_out/shapeit*", wk_dir))
        system(sprintf("rm -fr %s/mega/HapCut/phased*/*htrans_model", wk_dir))
        system(sprintf("rm -fr %s/mega/HapCut/hic_phased/generated_random_mix/*opt.*at.txt", wk_dir))
        system(sprintf("rm -fr %s/mega/aligned/merged_nodup_chrXchrY/*", wk_dir))

        # system(sprintf("rm -f %s/run*/hic_separated/*", wk_dir)) ## for hapcut
        # system(sprintf("rm -f %s/run*/splits/*n_sorted.REF*.bam", wk_dir)) ## for freebayes
        # system(sprintf("rm -f %s/run*/splits/hic_fixmate_sorted.bam*", wk_dir)) ## for hapcut
        system(sprintf('du -sh %s', wk_dir)) 

        if(mv2nas) 
        {
            out_log_mv_main = system(sprintf("mv %s %s", wk_dir, nas_dir), wait=TRUE)
            # out_log_mv_dump = system(sprintf("mv %s/%s* %s", scratch_dump_dir, LINE, nas_dump_dir), wait=TRUE)
            cat('out_log_mv_main: ', out_log_mv_main, '\n')
            # cat('out_log_mv_dump: ', out_log_mv_dump, '\n')
        }
    }

   get_hic_mat2plot = function(LINE, types, chr_num2look, range_start, range_end, bin_size)
    {
        get_cMat = function(mat, condition, index2keep=NULL)
        {
            mat = as.matrix(mat)
            density = apply(mat, 1, function(v) mean(v!=0))
            index_high = which(density > 0.05)

            if(is.null(index2keep)) index4cor = index_high
            if(!is.null(index2keep)) index4cor = intersect(index_high, index2keep)

            cMat_tmp = fast_cor(mat[index4cor, index4cor])
            ccMat_tmp = fast_cor(cMat_tmp)
            cMat = ccMat = mat*0
            cMat[index4cor, index4cor] = cMat_tmp
            ccMat[index4cor, index4cor] = ccMat_tmp
            # diag(cMat) = diag(ccMat) = 0
            li = list(cMat=cMat, ccMat=ccMat)
            return(li[[condition]])
        }

        get_hic_mat2plot_helper = function(LINE, mate, impute=TRUE, bin_size, norm=c('NONE', 'KR')[2], condition=c('observed', 'oe', 'cMat', 'ccMat')[2], chr_num2look, range_start, range_end)
        {
            if(mate=='bulk') hic_f = get_bulk_hic(LINE)
            if(mate %in% mates) hic_f = get_mate_hic(LINE, mate, impute=impute)
            ob_oe = condition
            if(condition %in% c('cMat', 'ccMat')) ob_oe = 'oe'

            bin_start = floor(range_start / bin_size) + 1
            bin_end = ceiling(range_end / bin_size)
            bin_names = sprintf('%s-%s', (bin_start:bin_end-1)*bin_size + 1, (bin_start:bin_end)*bin_size)
            index2keep = bin_start:bin_end

            chr2strawr = as.character(chr_num2look)
            mat_dump_deleted = try(silent=TRUE, strawr::straw(matrix=ob_oe, binsize=bin_size, chr1loc=chr2strawr, chr2loc=chr2strawr, unit='BP', fname=hic_f, norm=norm) %>% data.table)
            colnames(mat_dump_deleted) = c("pos_1", "pos_2", "counts")
            mat_dump_deleted = mat_dump_deleted[!is.nan(counts)]
            mat = convert_contact_format_fun(mat_dump_deleted, from='format_dump', to='format_sparse', chr_num=chr_num2look, binsize=bin_size)$format_sparse

            if(condition %in% c('cMat', 'ccMat')) 
            {
                corMat = get_cMat(mat, condition, index2keep=index2keep)
            }

            mat2plot = as.matrix(mat[bin_start:bin_end, bin_start:bin_end])
            if(condition %in% c('cMat', 'ccMat')) mat2plot = corMat[index2keep, index2keep]
            rownames(mat2plot) = colnames(mat2plot) = bin_names
            return(mat2plot)
        }

        mat2plot_li = list()
        mat2plot_li$parameters = data.table(ID=LINE, map_type=paste0(types, collapse=','), chr_num=chr_num2look, range_start=range_start, range_end=range_end, normalization='observed, oe (observed/expected), cor(oe), cor(cor(oe))')
        conditions = c('observed', 'oe', 'cMat', 'ccMat')

        mat_li = foreach(type=types) %do%
        {
            mat_li = foreach(condition=conditions) %do% get_hic_mat2plot_helper(LINE, mate=type, bin_size=bin_size, condition=condition, chr_num2look=chr_num2look, range_start=range_start, range_end=range_end)
            names(mat_li) = conditions
            mat_li
        }

        names(mat_li) = types
        mat2plot_li = c(mat2plot_li, mat_li)
        return(mat2plot_li)
    }      

    ################################################# add to mega for given runs

    combine_merge_nodup_runs = function(wk_dir, run_nums2combine)
    {
        chr_pairs = rbind(chr_pairs_intra, chr_pairs_inter)
        run_names2combine = paste0("run", run_nums2combine)
        split_dir_mega = paste0(file.path(wk_dir, 'mega', 'aligned', 'merged_nodup_chrXchrY'), "_run", paste0(run_nums2combine, collapse="_"))
        dir.create(split_dir_mega, recursive=TRUE, showWarnings=FALSE)

        no_out = foreach(i=1:nrow(chr_pairs)) %do%
        {
            chr_num_x = chr_pairs[[1]][i]
            chr_num_y = chr_pairs[[2]][i]

            split_dir_runs = file.path(wk_dir, run_names2combine, 'aligned', 'merged_nodup_chrXchrY')
            merged_nodup_runs = paste0(sprintf('%s/merged_nodups.txt_chr%schr%s.snps', split_dir_runs, chr_num_x, chr_num_y), collapse=' ')
            merged_nodup_mega = sprintf('%s/merged_nodups.txt_chr%schr%s.snps', split_dir_mega, chr_num_x, chr_num_y)

            system(sprintf("cat %s > %s", merged_nodup_runs, merged_nodup_mega))

            n_lines_runs = as.numeric( strsplit(trimws( tail(system( sprintf("wc -l %s", merged_nodup_runs), intern=TRUE ), 1) ), ' ')[[1]][1] ) ## total lines
            n_lines_mega = as.numeric( strsplit(trimws( tail(system( sprintf("wc -l %s", merged_nodup_mega), intern=TRUE ), 1) ), ' ')[[1]][1] ) ## total lines
            if(n_lines_runs != n_lines_mega) stop('n_lines_runs != n_lines_mega')
            cat(i, ' done\n')
        }
    }
