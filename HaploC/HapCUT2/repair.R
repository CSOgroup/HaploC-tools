    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    wk_dir = args[1]
    chunk_num = args[2]

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
    source(sprintf('%s/header/HaploC_header.R', repo_dir))

    #########################################################

    if(genome_version=='mm10') chr_nums = c(1:19, 'X')
    if(genome_version=='hg19') chr_nums = c(1:22, 'X')
    chr_chars = paste0('chr', chr_nums)
    n_chr = length(chr_nums)

    #########################################################

    wk_dir_chunk    = sprintf('%s/chunk%s/', wk_dir, chunk_num)  
    tmp_dir       = paste0(wk_dir_chunk, 'hapcut_tmp/')
    hic_split_dir = paste0(wk_dir_chunk, 'hic_separated/')
    fastq_prefix  = gsub('_R2.n_sorted.bam', '', gsub('_R1.n_sorted.bam', '', list.files(file.path(wk_dir_chunk, '/splits'), 'n_sorted.bam')))[1]
   
    #########################################################

    dir.create(hic_split_dir, showWarnings=FALSE)
    dir.create(tmp_dir, showWarnings=FALSE)
    
    #########################################################

    input_mate1 = paste0(wk_dir_chunk, '/splits/', fastq_prefix, '_R1.n_sorted.bam')
    input_mate2 = paste0(wk_dir_chunk, '/splits/', fastq_prefix, '_R2.n_sorted.bam')
       
    hic_repaired = paste0(wk_dir_chunk, '/splits/hic_repaired.sam')
    hic_fixmate_sorted = paste0(wk_dir_chunk, '/splits/hic_fixmate_sorted.bam')

    ######################################################### repair chimeric, use conda activate HapCUT2

    # repair chimeric Hi-C read pairs (from off-center Hi-C linker)
    command_repair = sprintf('python2.7 %s %s %s %s %s', HiC_repair_py, input_mate1, input_mate2, hic_repaired, file.path(repo_dir, "HaploC/dependancies/HapCUT2/"))
    system(command_repair)

    ######################################################### fixmate, use conda activate HaploC
 
    command_fixmate = sprintf('%s fixmate -@ 20 %s - | %s sort -@ 20 -T %s -o %s - && rm %s', samtools, hic_repaired, samtools, tmp_dir, hic_fixmate_sorted, hic_repaired)
    command_index =  sprintf('samtools index -@ 20 %s', hic_fixmate_sorted)
    system(command_fixmate)
    system(command_index)

    ######################################################### split into chrs

    no_out = foreach(chr_char=chr_chars) %do%
    {
        command_subset = sprintf('%s view -b %s %s -@ 10 -o %s/hic_%s.bam', samtools, hic_fixmate_sorted, chr_char, hic_split_dir, chr_char) 
        system(command_subset)
        cat(chr_char, ' split done\n')
    }

    #########################################################

    # fz = file.size( sprintf('%s/hic_%s.bam', hic_split_dir, chr_char) )/1024/1024 ## Mb
    # # if(fz > 1)
    # {
    #     unlink(input_mate1)
    #     unlink(input_mate2) 
    # }