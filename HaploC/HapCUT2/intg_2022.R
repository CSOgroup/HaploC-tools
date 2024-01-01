    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    wk_dir = args[1]
    chunks2consider = args[2] ## format: 1,2,4,5,7
    chr_num = args[3]

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

    #########################################################

    source(sprintf('%s/header/HaploC_header.R', repo_dir))
    registerDoParallel(cores=7) ## for one chr only
   
    ################################################################################

    chr_char = paste0('chr', chr_num)
    chunks = strsplit(chunks2consider, ",")[[1]]
    chunk_names = paste0('chunk', chunks)
    VCF_DIR      = paste0(wk_dir, '/mega/VCFs/')
    VCF_prefix   = 'mega'
    ref_fa = ref_fa_li[[genome_version]]

    dir.create(paste0(wk_dir, '/mega/HapCut/fragment_hapcut'), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(wk_dir, '/mega/HapCut/phased_hap_hapcut'), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(wk_dir, '/mega/HapCut/hic_phased'), recursive=TRUE, showWarnings=FALSE)      
    if(genome_version=='hg19') dir.create(paste0(wk_dir, '/mega/HapCut/phased_hap_intg'), recursive=TRUE, showWarnings=FALSE)
    if(genome_version=='hg19') dir.create(paste0(wk_dir, '/mega/HapCut/fragment_intg'), recursive=TRUE, showWarnings=FALSE)
    if(genome_version=='hg19') dir.create(paste0(wk_dir, '/mega/HapCut/shapeit_out'), recursive=TRUE, showWarnings=FALSE) 

    ################################################################################

    # runs = strsplit(run_nums2consider, "_")[[1]]
    runs = 1:2
    run_names = paste0('chunk', runs)
    VCF_DIR      = paste0(wk_dir, '/mega/VCFs/')
    VCF_prefix   = 'mega'
    
    dir.create(paste0(wk_dir, '/mega/HapCut/fragment_hapcut'), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(wk_dir, '/mega/HapCut/fragment_intg'), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(wk_dir, '/mega/HapCut/phased_hap_hapcut'), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(wk_dir, '/mega/HapCut/phased_hap_intg'), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(wk_dir, '/mega/HapCut/phased_contacts'), recursive=TRUE, showWarnings=FALSE)      
    dir.create(paste0(wk_dir, '/mega/HapCut/shapeit_out'), recursive=TRUE, showWarnings=FALSE) 

    ################################################################################

    shapeit2_fun = function(chr_char, QUAL_thresh)
    {
        frag_hapcut_chr = 'not_use'
        frag_intg_chr = 'not_use'

        vcf_chr_raw = sprintf('%s/%s.%s.vcf', VCF_DIR, VCF_prefix, chr_char)
        vcf_chr = gsub("HET", "ALL", vcf_filter(vcf_chr_raw, QUAL_thresh=QUAL_thresh, name_only=TRUE)) ## keep both HET and HOMO, i.e., this returns: mega.chr1.ALL.bi.QUAL=1.vcf

        hap_file = paste0(file.path(wk_dir, '/mega/HapCut/shapeit_out/', basename(vcf_chr)), '.haps')
        # ################################### shapeit2
        command_shapeit2 = sprintf('%s %s %s %s %s %s %s %s', SHAPEIT2, gsub('chr', '', chr_char), vcf_chr, frag_hapcut_chr, frag_intg_chr, paste0(wk_dir, '/mega/HapCut/shapeit_out'), 'create', repo_dir) 
        system(command_shapeit2)
        vcf2tab(hap_file, format_shapeit2_hap=TRUE) ## save the phased hap into wig
    }

    intg_fun_mini = function(run_names, chr_char, QUAL_thresh, maxIS_Mb, skip_if_exist=TRUE)
    {
        frag_hapcut_chr = sprintf('%s/mega/HapCut/fragment_hapcut/hic_%s_QUAL=%s_maxIS=%s.fragment', wk_dir, chr_char, QUAL_thresh, maxIS_Mb)
        frag_intg_chr = sprintf('%s/mega/HapCut/fragment_intg/hic_%s_QUAL=%s_maxIS=%s.fragment', wk_dir, chr_char, QUAL_thresh, maxIS_Mb)

        vcf_chr_raw = sprintf('%s/%s.%s.vcf', VCF_DIR, VCF_prefix, chr_char)
        vcf_chr = vcf_filter(vcf_chr_raw, QUAL_thresh=QUAL_thresh, name_only=TRUE) ## this line should not be included to the full dopar loop, with be overwritten
        vcf_chr = gsub("HET", "ALL", vcf_chr) ## because using ALL instead of HET. Otherwise the program will break I think

        extracthairs_runs = function(run_names) ## combine different runs together
        {
            for(run_num in 1:length(run_names))
            {
                run_name = run_names[run_num]
                bam_chr = sprintf('%s/%s/hic_separated/hic_%s.bam', wk_dir, run_name, chr_char) 

                ## --mbq is 13 by default -10*log10(0.05) = 13
                if(run_num==1) command_fragment = sprintf('%s --HiC 1 --bam %s --VCF %s --maxIS %s --indels 1 --mmq 10 --mbq 10 --ref %s > %s', EXTRACTHAIRS, bam_chr, vcf_chr, as.numeric(maxIS_Mb)*1E6, ref_fa, frag_hapcut_chr) 
                if(run_num >1) command_fragment = sprintf('%s --HiC 1 --bam %s --VCF %s --maxIS %s --indels 1 --mmq 10 --mbq 10 --ref %s >> %s', EXTRACTHAIRS, bam_chr, vcf_chr, as.numeric(maxIS_Mb)*1E6, ref_fa, frag_hapcut_chr) 
                system(command_fragment)
            }
        }

        extracthairs_runs(run_names)

        output_hap_intg = sprintf('%s/mega/HapCut/phased_hap_intg/hic_%s_QUAL=%s_maxIS=%s.hap', wk_dir, chr_char, QUAL_thresh, maxIS_Mb)
        fz = file.size(output_hap_intg)/ 1E6 ## Mb
        if(!is.na(fz) & (fz > 0.1) & skip_if_exist)
        {
            cat('File exist: ', output_hap_intg, '\n')
            return(NA) ## IF file is aleady there, do not redo by default
        }

        {
            # ################################### hapcut
            output_hap = sprintf('%s/mega/HapCut/phased_hap_hapcut/hic_%s_QUAL=%s_maxIS=%s.hap', wk_dir, chr_char, QUAL_thresh, maxIS_Mb)  
            output_model = sprintf('%s/mega/HapCut/phased_hap_hapcut/hic_%s_QUAL=%s_maxIS=%s.htrans_model', wk_dir, chr_char, QUAL_thresh, maxIS_Mb)  
            command_hapcut = sprintf('%s --fragments %s --vcf %s --output %s --hic 1 --htrans_data_outfile %s --outvcf 1', HapCUT22, frag_hapcut_chr, vcf_chr, output_hap, output_model)
            print(command_hapcut)
            system(command_hapcut)

            # ################################### shapeit2 append psudo reads
            command_shapeit2 = sprintf('%s %s %s %s %s %s %s', SHAPEIT2, gsub('chr', '', chr_char), vcf_chr, frag_hapcut_chr, frag_intg_chr, paste0(wk_dir, '/mega/HapCut/shapeit_out'), 'append') 
            system(command_shapeit2)

            # ################################### shapeit2 + hapcut = intg
            output_model = sprintf('%s/mega/HapCut/phased_hap_intg/hic_%s_QUAL=%s_maxIS=%s.htrans_model', wk_dir, chr_char, QUAL_thresh, maxIS_Mb) 
            command_intg = sprintf('%s --fragments %s --vcf %s --output %s --hic 1 --htrans_data_outfile %s --outvcf 1', HapCUT22, frag_intg_chr, vcf_chr, output_hap_intg, output_model)
            system(command_intg)
        }
    }

    logfile = paste0(wk_dir, sprintf('/mega/intg_log.chr%s.txt', chr_num))
    system(sprintf('date > %s', logfile))

    ################################################################################
    
    QUAL_threshs = 10
    maxIS_Mbs = c("50", '100') ## I decide to remove 1000, as it almost never appear in the opt set

    chr_chars = paste0('chr', chr_num)

    para_li = list(chr_char=chr_chars, maxIS_Mb=maxIS_Mbs, QUAL_thresh=QUAL_threshs)
    para_tab = data.table(do.call(expand.grid, para_li))
    para_tab$chr_char = as.vector(para_tab$chr_char)
    for(i in 1:ncol(para_tab)) para_tab[[i]] = as.vector(para_tab[[i]])

    para_tab_shapeit = unique(para_tab[, mget(c('chr_char', 'QUAL_thresh'))])

    ################################################################################ shapeit

    # no_out = foreach(k=1:nrow(para_tab_shapeit)) %dopar%
    # {
    #     chr_char = para_tab_shapeit$chr_char[k]
    #     QUAL_thresh = para_tab_shapeit$QUAL_thresh[k]
    #     try(shapeit2_fun(chr_char, QUAL_thresh))      
    # }


    ################################################################################ intg

    no_out = foreach(k=1:nrow(para_tab)) %dopar%
    {
        chr_char = para_tab$chr_char[k]
        QUAL_thresh = para_tab$QUAL_thresh[k]
        maxIS_Mb = para_tab$maxIS_Mb[k]

        try(intg_fun_mini(run_names=run_names, chr_char, QUAL_thresh, maxIS_Mb, skip_if_exist=FALSE))
        cat(sprintf('n_total=%s | k=%s, %s, QUAL_thresh:%s, maxIS_Mb:%s | done at: ', nrow(para_tab), k, chr_char, QUAL_thresh, maxIS_Mb), file=logfile, append=TRUE)
        system(sprintf('date >> %s', logfile))        
    }

    ################################################################################ phasing stat









    # phased_plot_dir = sprintf("%s/mega/plots/", wk_dir)
    # dir.create(phased_plot_dir, recursive=TRUE, showWarnings=FALSE)

    # para_li = list(maxIS_Mb=maxIS_Mbs, chr_char=chr_char, QUAL_thresh=QUAL_threshs)
    # para_tab_plot = data.table(do.call(expand.grid, para_li))
    # for(i in 1:ncol(para_tab_plot)) para_tab_plot[[i]] = as.vector(para_tab_plot[[i]])

    # ps = list()

    # for(k in 1:nrow(para_tab_plot))
    # {
    #     chr_char = para_tab_plot$chr_char[k]
    #     QUAL_thresh = para_tab_plot$QUAL_thresh[k]
    #     maxIS_Mb = para_tab_plot$maxIS_Mb[k]

    #     phased_vcf = sprintf("%s/mega/HapCut/phased_hap_intg/hic_%s_QUAL=%s_maxIS=%s.hap.phased.VCF", wk_dir, chr_char, QUAL_thresh, maxIS_Mb)
    #     ## Above should fileter out the homo SNPs to do the computing

    #     command_ps_stat = sprintf('bcftools filter %s -i \'GT!="0/0" && GT!="1/1"\' | bcftools query -f "[%%PS]\n" | sort | uniq -c', phased_vcf) ## phase set
    #     ps_stat_raw = trimws(system(command_ps_stat, intern=TRUE))
    #     ps_stat = data.table(do.call(rbind, strsplit(ps_stat_raw, " ")))
    #     colnames(ps_stat) = c("count", "PS")

    #     ps_stat$count = as.numeric(ps_stat$count)

    #     ps_stat = rbind(ps_stat[PS=="."], ps_stat[PS!="."][order(-count)])
    #     ps_stat$ratio = ps_stat$count / sum(ps_stat$count)
    #     ps_stat$index = 1:nrow(ps_stat)
    #     ps_stat = ps_stat[1:21, ]

    #     p = ggplot(data=ps_stat, aes(x=index, y=ratio, label=substr(ratio, 1, 5))) + geom_line(color="red") + geom_point()
    #     p = p + geom_label() + scale_x_continuous("chr", labels=c("!phased", 1:20), breaks=1:21) + labs(title=sprintf("%s:QUAL=%s:IS=%s", chr_char, QUAL_thresh, maxIS_Mb))
    #     ps[[k]] = p + ylim(0,1)
    # }

    # pdf(sprintf("%s/%s", phased_plot_dir, "phased_hap_intg_ratio.pdf"), width=8, height=6)
    # for(p in ps) print(p)
    # dev.off()

    # ############################

    # ps = list()

    # for(k in 1:nrow(para_tab_plot))
    # {
    #     chr_char = para_tab_plot$chr_char[k]
    #     QUAL_thresh = para_tab_plot$QUAL_thresh[k]
    #     maxIS_Mb = para_tab_plot$maxIS_Mb[k]

    #     phased_vcf = sprintf("%s/mega/HapCut/phased_hap_hapcut/hic_%s_QUAL=%s_maxIS=%s.hap.phased.VCF", wk_dir, chr_char, QUAL_thresh, maxIS_Mb)
    #     command_ps_stat = sprintf('bcftools filter %s -i \'GT!="0/0" && GT!="1/1"\' | bcftools query -f "[%%PS]\n" | sort | uniq -c', phased_vcf) ## phase set
    #     ps_stat_raw = trimws(system(command_ps_stat, intern=TRUE))
    #     ps_stat = data.table(do.call(rbind, strsplit(ps_stat_raw, " ")))
    #     colnames(ps_stat) = c("count", "PS")

    #     ps_stat$count = as.numeric(ps_stat$count)

    #     ps_stat = rbind(ps_stat[PS=="."], ps_stat[PS!="."][order(-count)])
    #     ps_stat$ratio = ps_stat$count / sum(ps_stat$count)
    #     ps_stat$index = 1:nrow(ps_stat)
    #     ps_stat = ps_stat[1:21, ]

    #     p = ggplot(data=ps_stat, aes(x=index, y=ratio, label=substr(ratio, 1, 5))) + geom_line(color="red") + geom_point()
    #     p = p + geom_label() + scale_x_continuous("chr", labels=c("!phased", 1:20), breaks=1:21) + labs(title=sprintf("%s:QUAL=%s:IS=%s", chr_char, QUAL_thresh, maxIS_Mb))
    #     ps[[k]] = p + ylim(0,1)
    # }

    # pdf(sprintf("%s/%s", phased_plot_dir, "phased_hap_hapcut_ratio.pdf"), width=8, height=6)
    # for(p in ps) print(p)
    # dev.off()


    # ############################

    # foreach(chr_num=chr_nums) %do%
    # {
    #     command = sprintf("Rscript %s/HaploC/stat_shapeit_hapcut_intg.R %s %s", repo_dir, wk_dir, chr_num)
    #     system(command)
    # }













