    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    wk_dir = args[1]

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

	######################################################### choose the number of chunks to get at least 40x depth

	get_chunks2use = function(x, n_cores)
	{
		x = as.numeric(x) ## 30X
		n_reads = foreach(chunk=1:n_chunks, .combine=c) %dopar%
		{
			splits_dir = sprintf('%s/chunk%s/splits', wk_dir, chunk)
			R1_bam = list.files(splits_dir, recursive=TRUE, full.names=TRUE, 'R1.n_sorted.bam')
			cmd_count = sprintf('samtools view -c %s', R1_bam)
			n_reads = as.numeric(system(cmd_count, intern=TRUE))
			n_reads
		}

		n_reads_tab = data.table(n=n_reads, chunk=1:n_chunks, n_total=round(cumsum(n_reads)/1E6, 1))
		chunk2look = which.max(n_reads)

		splits_dir = sprintf('%s/chunk%s/splits', wk_dir, chunk2look)
		R1_bam = list.files(splits_dir, recursive=TRUE, full.names=TRUE, 'R1.n_sorted.bam')
		R1_cov = gsub('n_sorted.bam', 'cov', R1_bam)
		system(sprintf("samtools sort -T %s -@ %s %s -O BAM | samtools coverage - > %s ", splits_dir, n_cores, R1_bam, R1_cov))
		cov_chrs = fread(R1_cov)
		colnames(cov_chrs)[1] = 'chr'
		med_cov = median(cov_chrs[chr %in% chr_chars[1:19]]$meandepth) ## keep only the first 19 chrs will be enough

		n_reads_tab$cov = 0
		n_reads_tab[chunk2look]$cov = med_cov
		min_reads = x / med_cov * n_reads_tab$n[chunk2look] / 1E6 ## set the min coverage as 40X
		if(max(n_reads_tab$n_total) < min_reads) min_chunk = n_chunks
		if(max(n_reads_tab$n_total) >= min_reads) min_chunk = min(n_reads_tab[n_total > min_reads]$chunk)

		cat('n chunks to use:', min_chunk, '\n')

		chunks = 1:min_chunk
		return(chunks)
	}

	#########################################################

    ref_fa = ref_fa_li[[genome_version]]

    if(genome_version=='mm10') chr_nums = c(1:19, 'X')
    if(genome_version=='hg19') chr_nums = c(1:22, 'X')
    chr_chars = paste0('chr', chr_nums)
    n_chr = length(chr_nums)
        
    ######################################################### cannot use samtools view here, because n_sorted_bam is not sorted
   
    stopImplicitCluster()
    n_cores = 23
    registerDoParallel(cores=n_cores)

    #########################################################
    
    n_chunks = fread(sprintf('%s/log_file/n_chunks.txt', wk_dir))[[1]]
	chunks = get_chunks2use(x, n_cores)

    ######################################################### split by chr

    n_sorted_bams = sort(list.files(file.path(wk_dir, paste0("chunk", chunks)), recursive=TRUE, full.names=TRUE, 'n_sorted.bam'))
    no_out = foreach(n_sorted_bam=n_sorted_bams) %dopar% ## split into individual chrs
    {
        cmd = sprintf("%s split -in %s -reference", bamtools, n_sorted_bam)
        cat(cmd, "\n")
        system(cmd, wait=TRUE)
    }

    #########################################################

    splits_dir = file.path(wk_dir, 'mega/splits/')
    vcf_dir = file.path(wk_dir, 'mega/VCFs/')
    dir.create(recursive=TRUE, splits_dir, showWarnings=FALSE)
    dir.create(recursive=TRUE, vcf_dir, showWarnings=FALSE)

    #########################################################

    no_out = foreach(chr_num=rev(chr_nums)) %dopar% ## do for each chr. It is fine to do it in parallel, will not take too much memory. Otherwise it takes too long / change back to not parallel
    {
        mega_bam = file.path(wk_dir, 'mega/splits/', sprintf('mega.pos_sorted.chr%s.bam', chr_num))
        n_sorted_bams_chr = paste0(gsub("n_sorted.bam", sprintf("n_sorted.REF_chr%s.bam", chr_num), n_sorted_bams), collapse=" ")
        # cmd_merge = sprintf('%s merge -f -@ n_chr -n %s %s', samtools, paste0(mega_bam, '.tmp'), n_sorted_bams_chr)
        # cat(cmd_merge, file=gsub("bam", "log", mega_bam))

        # cmd_merge_sort = sprintf('%s sort -@ n_chr -m 2G %s -O BAM -o %s -T %s/tmp_chr%s', samtools, paste0(mega_bam, '.tmp'), mega_bam, splits_dir, chr_num)
        cmd_merge_sort = sprintf('samtools merge -n -f - %s | samtools sort - -O BAM -o %s -T %s/tmp_chr%s', n_sorted_bams_chr, mega_bam, splits_dir, chr_num)            
        cmd_index =  sprintf('%s index %s', samtools, mega_bam)

        cat(cmd_merge_sort, file=gsub("bam", "log", mega_bam))
        system(cmd_merge_sort)

        # system(cmd_merge)
        system(sprintf("rm %s", n_sorted_bams_chr))            
        # unlink(paste0(mega_bam, '.tmp'))
        system(cmd_index)
    }

    ######################################################### call variants

    bam_head = file.path(wk_dir, 'mega/splits/', sprintf('mega.pos_sorted'))
    parallel_free_v2(free_exe, ref_fa, bam_head=bam_head, wk_dir=wk_dir, chrs=chr_nums) ## v0.9.21 before Nov 2023
    
    ######################################################### stat

    n_snp_file = file.path(wk_dir, 'mega/VCFs/', 'n_variants_het.no_ins.txt')
    n_snp_plot_file = file.path(wk_dir, 'mega/VCFs/', 'n_variants_het.no_ins.pdf')
    QUAL_plot_file = file.path(wk_dir, 'mega/VCFs/', 'QUAL_distr_het.no_ins.pdf')

    vcf_prefix = gsub('pos_sorted', '', bam_head)

    no_out = foreach(chr=chr_chars, .combine=rbind) %do% ## add more QS choices, 21 Dec 2022
    {
        vcf_des_file = gsub('splits', 'VCFs', paste0(vcf_prefix, chr, '.vcf'))
        vcf_filter(vcf_raw_file=vcf_des_file, QUAL_thresh=10, name_only=FALSE, heto_only=FALSE) ## add more QS choices
        vcf_filter(vcf_raw_file=vcf_des_file, QUAL_thresh=20, name_only=FALSE, heto_only=FALSE)      
    }

    n_snps = foreach(QUAL_thresh=c(10,20), .combine=rbind) %do%
    {
        n_snps = foreach(chr=chr_chars, .combine=rbind) %do% ## pthread_create failed for thread 68 of 88: Resource temporarily unavailable
        {
            vcf_des_file = gsub('splits', 'VCFs', paste0(vcf_prefix, chr, sprintf('.HET.bi.QUAL=%s.vcf', QUAL_thresh)))
            count = system( sprintf('%s view -H %s | wc -l', bcftools, vcf_des_file), intern=TRUE )
            c(chr, count)
        }
        tab = data.table(n_snps, QUAL=QUAL_thresh)
        tab[[1]] = 1:length(chr_nums)
        tab
    }

    n_snps[[2]] = as.numeric(n_snps[[2]])
    colnames(n_snps)[1:2] = c('chr', 'count')
    
    #########################################################

    fwrite(n_snps, file=n_snp_file, row.names=FALSE, col.names=FALSE, sep='\t')  
   
    #########################################################

    ps = foreach(QUAL_thresh=10) %do%
    {
        df = n_snps[QUAL==QUAL_thresh]
        p = ggplot(data=df, aes(x=chr, y=count, label=count)) + geom_line(color="red") + geom_point()
        p = p + geom_label() + scale_x_continuous("chr", labels=chr_nums, breaks=1:n_chr)
        p + coord_cartesian(ylim=c(0, max(n_snps$count)*1.1))
    }

    p = ggplot(data=n_snps, aes(x=chr, y=count)) + geom_line(aes(linetype=as.factor(QUAL), color=as.factor(QUAL))) + geom_point()
    p = p + scale_x_continuous("chr", labels=chr_nums, breaks=1:n_chr)
    p = p + coord_cartesian(ylim=c(0, max(n_snps$count)*1.1)) + theme(legend.position='none')

    ps[[2]] = p

    pdf(n_snp_plot_file, width=8, height=6)
    for(p in ps) print(p)    
    graphics.off()

    #########################################################

	pdf(QUAL_plot_file)
	for(chr in chr_chars)
	{
	    vcf_des_file = gsub('splits', 'VCFs', paste0(vcf_prefix, chr, '.HET.bi.QUAL=10.vcf'))
	    QUALs = sort(as.numeric(system(sprintf('%s query -f "%%QUAL\n" %s', bcftools, vcf_des_file), intern=TRUE)))
	    if(length(QUALs) < 10) 
	    {
	        plot(NULL, xlim=c(1,1000), col='blue', xlab='QUAL', ylab='ecdf', main=chr, xaxt = "n", ylim=c(0,0.8))
	        next
	    }

	    x = ecdf((QUALs))
	    plot(x, xlim=c(1,1000), col='blue', xlab='QUAL', ylab='ecdf', main=chr, xaxt = "n", ylim=c(0,0.8))
	    axis(1, at=c(10,100*(1:10)), labels=c(10,100*(1:10)), cex.axis=0.7)
	    abline(v=c(1,10,100), col='red', lty=2)
	}
	dev.off()

	#########################################################

	get_chr2use(wk_dir, genome_version=genome_version)    

	#########################################################