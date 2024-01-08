    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    wk_dir = args[1]
    job_name = args[2]

    ########################################################

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

    ########################################################

    if(length(enzyme)==0) {cat('enzyme not found\n'); quit(save='no')}

    ######################################################## -1: organize data

    if(job_name=='split')
    {
        data_num = as.numeric(args[3])
        
        sizeGb = as.numeric(sizeGb)
        data_source_f = sprintf('%s/log_file/data_source.txt', wk_dir)
        fastq_id = fread(data_source_f, header=FALSE)[[1]][data_num]
        cat('data:', fastq_id, '\n')
        if(is.na(sizeGb)) split_fastq(wk_dir, fastq_id, sizeGb=5)
        if(!is.na(sizeGb)) split_fastq(wk_dir, fastq_id, sizeGb=as.numeric(sizeGb))

        cat(sprintf('Finished: job=%s', job_name))
        quit(save="no")
    }

    if(job_name=='distr')
    {
        setwd(wk_dir)
        log_chunk_files = sprintf('%s/log_file/log.chunk_source.txt', wk_dir)

        fs_raw = list.files(wk_dir, recursive=TRUE, full.names=TRUE)
        fs = fs_raw[grepl('fastq.gz', fs_raw) & !grepl('fastq_bp', fs_raw)]
        fs = gsub('001.fastq.gz', '.fastq.gz', fs)
        fs = sort(fs)
        fs_bases = basename(fs)
        bases = substr(fs_bases, 1, nchar(fs_bases) - 11)
        if(substr(bases[1], nchar(bases[1]) - 1, nchar(bases[1]) - 1)=='_') bases = substr(fs_bases, 1, nchar(fs_bases) - 12)

        bases_uniq = unique(bases)
        n_chunks = length(bases_uniq)

        cat(paste0(sprintf('chunk%s:%s', 1:n_chunks, bases_uniq), collapse=' | '), file=log_chunk_files)

        for( i in 1:length(bases_uniq) )
        {
            system(sprintf('mkdir -p %s/chunk%s/fastq/', wk_dir, i))
            cmd_mv_R1 = paste0('mv ', bases_uniq[i], '_1.fastq.gz', ' chunk', i, '/fastq/', bases_uniq[i], '_R1.fastq.gz')
            cmd_mv_R2 = paste0('mv ', bases_uniq[i], '_2.fastq.gz', ' chunk', i, '/fastq/', bases_uniq[i], '_R2.fastq.gz')
            system(cmd_mv_R1)
            system(cmd_mv_R2)  
        }
        cat(n_chunks, file=sprintf('%s/log_file/n_chunks.txt', wk_dir))
        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }

    ######################################################## do clumpy: 0 / bwa: 1

    if(job_name %in% c('bwa', 'hic'))
    {
        chunk_num = as.numeric(args[3])

        #########################################################

        res_sites  = sprintf("%s/genomicData/restriction_sites/%s_%s.txt", repo_dir, genome_version, enzyme)
        chrom_sizes = sprintf("%s/genomicData/references/%s.chrom.sizes", repo_dir, genome_version)
        ref_fa = ref_fa_li[[genome_version]]

        #########################################################

        chunk_name = paste0("chunk", chunk_num)
        fs = list.files(sprintf("%s/%s/fastq", wk_dir, chunk_name), recursive=TRUE)[1]
        fastq_prefix = unique(substr(fs, 1, nchar(fs)-12))
        log_global = sprintf('%s/%s/log_global.txt', wk_dir, chunk_name)

        if(job_name=='bwa') STEPs = c(0:1)
        if(job_name=='hic') STEPs = c(2:5)
 
        for(STEP in STEPs)
        {

            # if(STEP==1) system(sprintf('echo -n START STEP %s at `date` ::  > %s', STEP, log_global))
            # if(STEP!=1) system(sprintf('echo -n START STEP %s at `date` ::  >> %s', STEP, log_global))

            try(run_juicer(wk_dir, run_name=chunk_name, fastq_prefix=fastq_prefix, STEP=STEP, thread=thread4bwa, juicerDir=juicerDir, enzyme=enzyme, res_sites=res_sites, ref_fa=ref_fa, chrom_sizes=chrom_sizes, server_capacity_low=0))
            # system(sprintf('echo FINISH at `date` >> %s', log_global))
        }

        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }

    if(job_name=='mega')
    {
        res_sites  = sprintf("%s/restriction_sites/%s_%s.txt", juicerDir, genome_version, enzyme)
        run_juicer_mega(wk_dir, juicerDir, enzyme, res_sites)
        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }

    ######################################################### 2: call variants

    if(job_name=='snp')
    {
        ref_fa = ref_fa_li[[genome_version]]        
        n_chunks = fread(sprintf('%s/log_file/n_chunks.txt', wk_dir))[[1]]
        chunks2consider = paste0(1:n_chunks, collapse="_")
        cmd_free = sprintf('Rscript %s/SNP/snp_calling.R %s %s %s %s', HaploC_dir, wk_dir, 'bam_not_exsis', chunks2consider, ref_fa)
        system(cmd_free, ignore.stdout=FALSE, wait=TRUE)
        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }

    ######################################################### 3: repair

    if(job_name=='rep')
    {
        chunk_num = as.numeric(args[3])
        chunk_num = as.numeric(args[3])        
        cmd_repair = sprintf('%s %s/HapCUT2/repair.R %s %s', Rscript, HaploC_dir, wk_dir, chunk_num)
        try(system(cmd_repair, wait=TRUE))
        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
  
    }

    ######################################################### 4: intg & stat

    if(job_name=='intg')
    {
        n_chunks = fread(sprintf('%s/log_file/n_chunks.txt', wk_dir))[[1]]
        chr_num = args[4]
        
        chunks2consider = paste0(1:n_chunks, collapse=",")      
        cmd_intg = sprintf('Rscript %s/HapCUT2/intg.R %s %s %s', HaploC_dir, wk_dir, chunks2consider, chr_num)
        cmd_stat = sprintf('Rscript %s/HapCUT2/phasing_consistency.R %s %s', HaploC_dir, wk_dir, chr_num)
        system(ignore.stdout=FALSE, cmd_intg, wait=TRUE)
        if(genome_version=='hg19') system(ignore.stdout=FALSE, cmd_stat, wait=TRUE)        
        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }


    # cmd_stat = sprintf('for chr_num in {1..22} X; do Rscript %s/HaploC/%s %s $chr_num & done', repo_dir, 'stat_shapeit_hapcut_intg.R', wk_dir)
    # system(cmd_stat, wait=FALSE)

    ######################################################### 5: intersect with SNPs & gzip merged_nodup

    if(job_name=='sec')
    {
        chunk_num = as.numeric(args[3])
        cmd_sec = sprintf('%s %s/SNP/intersect_with_snp.R %s %s', Rscript, HaploC_dir, wk_dir, chunk_num)
        try(system(cmd_sec, wait=TRUE))
        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }

    # add2mega = c(TRUE, FALSE)[1]
    # identity = c('.HET.bi.QUAL=10', '_ground_truth')[1]
    # source(file.path(repo_dir, 'HaploC', 'intg_split_merged_nodup.R'))
    # system(ignore.stdout=FALSE, sprintf('Rscript %s %s', file.path(repo_dir, 'HaploC', 'gzip_merged_nodup.R'), wk_dir), wait=FALSE)

    ######################################################### 6: opt

    if(job_name=='opt')
    {
        chr_nums2use = get_chr2use(wk_dir)
        n_chr = length(chr_nums2use)

        intra_chr_pairs = data.table(expand.grid(chr_nums2use, chr_nums2use))[Var1 == Var2]
        inter_chr_pairs = data.table(expand.grid(1:23, 1:23))[Var1 < Var2]
        inter_chr_pairs[inter_chr_pairs==23] = 'X'
        inter_chr_pairs = inter_chr_pairs[Var1 %in% chr_nums2use][Var2 %in% chr_nums2use]
        chr_pairs = rbind(intra_chr_pairs, inter_chr_pairs)

        n_chunks = fread(sprintf('%s/log_file/n_chunks.txt', wk_dir))[[1]]

        #########################

        split_dir_mega = file.path(wk_dir, 'mega', 'aligned', 'merged_nodup_chrXchrY')
        dir.create(split_dir_mega, recursive=TRUE, showWarnings=FALSE)

        add_chunk = function(chunk_num)
        {
            chunk_name = sprintf('chunk%s', chunk_num)
            no_out = foreach(i=1:nrow(chr_pairs)) %do%
            {
                chr_num_x = chr_pairs[[1]][i]
                chr_num_y = chr_pairs[[2]][i]

                split_dir_chunks = file.path(wk_dir, chunk_name, 'aligned', 'merged_nodup_chrXchrY')
                merged_nodup_chunks = paste0(sprintf('%s/merged_nodups.txt_chr%schr%s.snps', split_dir_chunks, chr_num_x, chr_num_y), collapse=' ')
                merged_nodup_mega = sprintf('%s/merged_nodups.txt_chr%schr%s.snps', split_dir_mega, chr_num_x, chr_num_y)
                if(chunk_num!=1) system(sprintf("cat %s >> %s", merged_nodup_chunks, merged_nodup_mega), wait=TRUE)
                if(chunk_num==1) system(sprintf("cat %s > %s", merged_nodup_chunks, merged_nodup_mega), wait=TRUE)                
            }
        }

        for(chunk_num in 1:n_chunks) try(add_chunk(chunk_num))

        #########################

        output_hap = file.path(wk_dir, '/mega/HapCut/phased_hap_hapcut/')
        output_hap_intg = file.path(wk_dir, '/mega/HapCut/phased_hap_intg/')

        if(genome_version=='mm10') ## shapeit2 not used
        {
            cmd_rm = sprintf('rm -rf %s', output_hap_intg); system(cmd_rm)
            cmd_cp = sprintf('cp -r %s %s', output_hap, output_hap_intg); system(cmd_cp)
        } 

        #########################

        registerDoParallel(cores=length(chr_nums2use))
        no_out = foreach(chr_num=chr_nums2use) %dopar%
        {
            cmd_opt = sprintf('Rscript %s/block_assemble/%s %s %s %s', HaploC_dir, 'block_assemble.R', wk_dir, chr_num, repo_dir)
            system(cmd_opt, wait=TRUE)
        }

        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }
  
    ######################################################### 8: diploid_map, intra / do not do parallel when using graviton

    if(job_name=='mates')
    {
        cmd_dip_impute = sprintf("%s %s %s %s %s", diploid_opt_intra_sh, wk_dir, "impute", repo_dir, genome_version)
        cmd_dip_n_impute = sprintf("%s %s %s %s %s", diploid_opt_intra_sh, wk_dir, "n_impute", repo_dir, genome_version)
        cmd_generate_random = sprintf("%s %s %s %s", generate_random_mixed_sh, wk_dir, repo_dir, genome_version)

        create_inter_map = function()
        {
            script_dip_inter = sprintf('%s/create_phased_maps/diploid_opt_inter_helper.sh', HaploC_dir)
            identity = c("ground_truth", "opt")[2]

            wk_dir_mega = file.path(wk_dir, 'mega')

            ################################### remomve the file first / have to

            unlink(sprintf("%s/HapCut/hic_phased/*inter*", wk_dir_mega))
            
            #################################################

            chr_nums2use = get_chr2use(wk_dir, genome_version)

            inter_chr_pairs = data.table(expand.grid(1:23, 1:23))[Var1 < Var2]
            if(genome_version=='hg19') inter_chr_pairs[inter_chr_pairs==23] = 'X'
            if(genome_version=='mm10') inter_chr_pairs[inter_chr_pairs==20] = 'X'            
            inter_chr_pairs = inter_chr_pairs[Var1 %in% chr_nums2use][Var2 %in% chr_nums2use]
           
            ###################################

            inter_chr_pairs = inter_chr_pairs[Var1 %in% chr_nums2use][Var2 %in% chr_nums2use]

            no_out = foreach(impute=c("impute", "n_impute")) %do% ## do not do it in parallel, complicated
            {
                no_out = foreach(i=1:nrow(inter_chr_pairs)) %do%
                {
                    chr_num_x = inter_chr_pairs[[1]][i]
                    chr_num_y = inter_chr_pairs[[2]][i]

                    command_xy = sprintf("%s %s %s %s %s %s %s %s %s", script_dip_inter, wk_dir_mega, chr_num_x, chr_num_y, identity, 'no', impute, repo_dir, genome_version)
                    system(command_xy)
                }

                ###################################

                command_create_hic = sprintf("%s %s %s %s %s %s %s %s %s", script_dip_inter, wk_dir_mega, "chrNA", "chrNA", identity, 'yes', impute, repo_dir, genome_version) ## chrNA: pseudo chr
                system(command_create_hic)
            }

            ###################################

            fs = list.files(sprintf("%s/HapCut/hic_phased/", wk_dir_mega), full.names=TRUE)
            fs = fs[!grepl('[.]hic', fs) & !grepl('ALL', fs)]
            unlink(fs)
        }

        system(cmd_dip_impute, wait=TRUE) ## this will not take too much time / let the cmd_dip_n_impute to wait
        system(cmd_dip_n_impute, wait=TRUE)
        system(cmd_generate_random, wait=TRUE)
        
        create_inter_map()
        cat(sprintf('Finished:: job=%s', job_name))
        quit(save="no")
    }
