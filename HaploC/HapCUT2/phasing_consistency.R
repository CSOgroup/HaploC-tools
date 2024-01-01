	# !/usr/bin/env Rscript
	args = commandArgs(trailingOnly=TRUE)
	wk_dir = args[1]
	chr_num = args[2]

	####################################

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

	####################################

    if(genome_version=='mm10') chr_nums = c(1:19, 'X')
    if(genome_version=='hg19') chr_nums = c(1:22, 'X')
    chr_chars = paste0('chr', chr_nums)
    n_chr = length(chr_nums)

	####################################

	QUAL_threshs = 10 ## do not inlcude 1 and 100 / the higher the qual_thresh, the sparser the fragments
    maxIS_Mbs = strsplit(maxIS, split=',')[[1]]

	phasing_stat_plot_file = sprintf('%s/mega/HapCut/phasing_stat/phasing_stat_chr%s.pdf', wk_dir, chr_num)
	dir.create(dirname(phasing_stat_plot_file), recursive=TRUE, showWarnings=FALSE)

	#################################### shapeit

	types = c("HET", "ALL")
	hap_ALL_HET = sprintf("%s/mega/HapCut/shapeit_out/mega.chr%s.%s.bi.QUAL=10.vcf.haps", wk_dir, chr_num, types)
	types = types[file.exists(hap_ALL_HET)]

	vcf_tab_li_shapeit = foreach(type=types) %do%
	{
		hap_file = sprintf("%s/mega/HapCut/shapeit_out/mega.chr%s.%s.bi.QUAL=10.vcf.haps", wk_dir, chr_num, type)
		vcf_tab = vcf2tab(hap_file, format_shapeit2_hap=TRUE)
		vcf_tab$chr = paste0('chr', vcf_tab$chr)
		vcf_tab[, c("chr", "POS", "hap_val")]
	}

	names(vcf_tab_li_shapeit) = paste0("shapeit_", types)

	#################################### intg

	vcf_files = sprintf('%s/mega/HapCut/phased_hap_intg/hic_chr%s_QUAL=10_maxIS=%s.hap.phased.VCF', wk_dir, chr_num, maxIS_Mbs) 
	maxIS_Mbs = maxIS_Mbs[file.exists(vcf_files)] ## only keep those exists (some parameter can fail)

	para_li = list(chr=chr_chars, QUAL_thresh=QUAL_threshs, maxIS_Mb=maxIS_Mbs)
	quality_tab = data.table(do.call(expand.grid, para_li))
	for(i in 1:ncol(quality_tab)) quality_tab[[i]] = as.vector(quality_tab[[i]])

	quality_tab$para = paste0(quality_tab$QUAL_thresh, '_', quality_tab$maxIS_Mb)

	####################################

	block_ratio_thresh = 0.01 ## it was 0.05 before 02.08.2021
	chr_char = paste0('chr', chr_num)
	quality_tab = quality_tab[chr==chr_char]
	n = nrow(quality_tab)

	#################################### compute consistency. Some parameters the overall selected block can still be low, does not matter

	vcf_tab_li_intg = foreach(i=1:n) %dopar%
	{
		QUAL_i = quality_tab$QUAL_thresh[i]
		max_IS_i = quality_tab$maxIS_Mb[i]
		vcf_file = sprintf('%s/mega/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.hap.phased.VCF', wk_dir, chr_num, QUAL_i, max_IS_i)
		vcf_tab = try(vcf2tab(block_ratio_thresh=block_ratio_thresh, vcf_file, which_blocks=NA, gz=FALSE, phased=TRUE, save2wig=TRUE, consider_ps=TRUE, wig_ground_truth=NA, format_shapeit2_hap=FALSE))
		# if(class(vcf_tab) == "try-error")  vcf_tab = data.table(1) ## in case when all block ratios are small
		vcf_tab[, c("chr", "POS", "hap_val")]
	}

	names(vcf_tab_li_intg) = paste0('p', quality_tab$para)

	vcf_het_input = sprintf('%s/mega/VCFs/mega.chr%s.HET.bi.QUAL=10.vcf', wk_dir, chr_num)
	n_total = as.numeric(system(sprintf('bcftools filter %s -i \'GT!="0/0" && GT!="1/1"\' | wc -l ', vcf_het_input), intern=TRUE))

	####################################

	vcf_tab_li = c(vcf_tab_li_intg, vcf_tab_li_shapeit)
	# for(i in 1:length(vcf_tab_li)) colnames(vcf_tab_li[[i]])[3] = names(vcf_tab_li)[i]

	relative_phased = sapply(vcf_tab_li, nrow) / max(sapply(vcf_tab_li, nrow))
	relative_phased['max_in_total'] = max(sapply(vcf_tab_li, nrow)) / n_total

	####################################

	n = length(vcf_tab_li)
	para_li = list(i=1:n, j=1:n)
	para_tab = data.table(do.call(expand.grid, para_li))
	para_tab = para_tab[i>j]	  

	mean_switch_mat = mean_identical_mat = mean_identical_mat_switched = matrix(0, n, n)

	identical_vals = foreach(k=1:nrow(para_tab), .combine=rbind) %dopar%
	{
		i = para_tab[[1]][k]
		j = para_tab[[2]][k]
		vcf_tab_x = vcf_tab_li[[i]][, mget(c('POS', 'hap_val'))]
		vcf_tab_y = vcf_tab_li[[j]][, mget(c('POS', 'hap_val'))]
		vcf_tab_xy = merge(vcf_tab_x, vcf_tab_y, by='POS')[order(POS)]
		identical = mean(sign(vcf_tab_xy[[3]])==sign(vcf_tab_xy[[2]]))
		c(i, j, max(identical, 1 - identical))
	}

	mean_identical_mat[identical_vals[,1:2]] = identical_vals[,3]
	mean_identical_mat = as.matrix(my_forceSymmetric(mean_identical_mat))


	switch_vals = foreach(k=1:nrow(para_tab), .combine=rbind) %dopar%
	{
		i = para_tab[[1]][k]
		j = para_tab[[2]][k]
		vcf_tab_x = vcf_tab_li[[i]][, mget(c('POS', 'hap_val'))]
		vcf_tab_y = vcf_tab_li[[j]][, mget(c('POS', 'hap_val'))]
		vcf_tab_xy = merge(vcf_tab_x, vcf_tab_y, by='POS')[order(POS)]
		vcf_tab_diff = abs(vcf_tab_xy[[3]] - vcf_tab_xy[[2]])/2
		c(i, j, mean(abs(diff(vcf_tab_diff)))) ## maybe to redefine this value using the official switch error
	}

	mean_switch_mat[switch_vals[,1:2]] = switch_vals[,3]
	mean_switch_mat = as.matrix(my_forceSymmetric(mean_switch_mat))

	rownames(mean_switch_mat) = colnames(mean_switch_mat) = rownames(mean_identical_mat) = colnames(mean_identical_mat) = names(vcf_tab_li)
	diag(mean_identical_mat) = 1

	####################################

	# phasing_stat = list(relative_phased=relative_phased, mean_identical_mat=mean_identical_mat, mean_switch_mat=mean_switch_mat)
	# R_f = file.path(wk_dir, '/mega/HapCut/phasing_stat/phasing_stat.Rds')
	# saveRDS(phasing_stat, file=R_f)

	#################################### make the plots

	pdf(phasing_stat_plot_file, width=9, height=9)
	plot(x=seq_along(relative_phased), relative_phased, type = "b", frame = FALSE, pch = 19, col = "red", ylim=c(0,1), main="relative_phased snps")
	box()
	abline(h=0.8, lty=3, col="blue")
	text(x=seq_along(relative_phased), y=relative_phased - 0.07, names(relative_phased), srt=90)

	cellnote = matrix(substr(mean_identical_mat, 1, 4), nrow(mean_identical_mat))
	x = (mean_identical_mat^(1/3))
	# rownames(x) = paste0(rownames(x), '( ', gsub('err', 'err', quality_tab$flag))
	heatmap.2(main="mean_identical_mat[0to1]", margins=c(12,12), revC=TRUE, key=0, cellnote=cellnote, x=x, dendrogram='none', Rowv=FALSE, Colv=FALSE, col="bluered", trace="none", ylab="", xlab="")


	cellnote = matrix(substr(mean_switch_mat, 1, 5), nrow(mean_switch_mat))
	x = (mean_switch_mat^(1/3))
	# rownames(x) = paste0(rownames(x), '( ', gsub('err', 'err', quality_tab$flag))
	heatmap.2(main="mean_switch_mat[0to1]", margins=c(12,12), revC=TRUE, key=0, cellnote=cellnote, x=x, dendrogram='none', Rowv=FALSE, Colv=FALSE, col="bluered", trace="none", ylab="", xlab="")
	dev.off()

	####################################
	
