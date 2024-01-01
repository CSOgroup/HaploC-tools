    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    wk_dir = args[1]
    chr_num = args[2]

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
    chr_nums2use = get_chr2use(wk_dir, genome_version) ## character
    if(!(chr_num %in% chr_nums2use))
    {
        # cat('SNP count low, do no phase\n', sprintf('Finished:: job=intg'))
        quit(save="no")
    } 
   
    ################################################################################

    QUAL_thresh = 10
    if(genome_version=='hg19' & chr_num=='23') chr_num='X'
    if(genome_version=='mm10' & chr_num=='20') chr_num='X'    
    chr_char = paste0('chr', chr_num)

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

    shapeit2_fun = function(chr_char, QUAL_thresh)
    {
        frag_hapcut_chr = 'not_use'
        frag_intg_chr = 'not_use'

        vcf_chr = sprintf("%s/mega.%s.ALL.bi.QUAL=%s.vcf", VCF_DIR, chr_char, QUAL_thresh)
        hap_file = paste0(file.path(wk_dir, '/mega/HapCut/shapeit_out/', basename(vcf_chr)), '.haps')
        
        #################################### shapeit2
        
        command_shapeit2 = sprintf('%s %s %s %s %s %s %s %s', SHAPEIT2, gsub('chr', '', chr_char), vcf_chr, frag_hapcut_chr, frag_intg_chr, paste0(wk_dir, '/mega/HapCut/shapeit_out'), 'create', repo_dir) 
        system(command_shapeit2)
        vcf2tab(hap_file, format_shapeit2_hap=TRUE) ## save the phased hap into wig
        return()
    }

    logfile = paste0(wk_dir, sprintf('/mega/intg_log.chr%s.txt', chr_num))
    system(sprintf('date > %s', logfile))

    ################################################################################ shapeit
    
    if(genome_version=='hg19') try(shapeit2_fun(chr_char, QUAL_thresh))

    ################################################################################ 