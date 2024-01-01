	# !/usr/bin/env Rscript

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
	source(sprintf('%s/header/HaploC_header.R', repo_dir))

	####################################

	dir.create(sprintf('%s/mega/HapCut/fragile_sites/', wk_dir), showWarnings=FALSE)

	fragile_flip_f = sprintf('%s/mega/HapCut/fragile_sites/fragile_sites_flip_chr%s.Rdata', wk_dir, chr_num)
	fragile_plot_f = sprintf('%s/mega/HapCut/fragile_sites/fragile_sites_chr%s.pdf', wk_dir, chr_num)
	consist_plot_f = sprintf('%s/mega/HapCut/fragile_sites/consist_heatmap_chr%s.pdf', wk_dir, chr_num)
	consist_rds_f = sprintf('%s/mega/HapCut/fragile_sites/consist_mat_chr%s.rds', wk_dir, chr_num)

	####################################

	QUAL_threshs = 10
	maxIS_Mbs = strsplit(maxIS, split=',')[[1]]

	vcf_files = sprintf('%s/mega/HapCut/phased_hap_intg/hic_chr%s_QUAL=10_maxIS=%s.hap.phased.VCF', wk_dir, chr_num, maxIS_Mbs)
	maxIS_Mbs = maxIS_Mbs[file.exists(vcf_files)] ## only keep those exists (some parameter can fail)

	chr_chars = paste0('chr', get_chr2use(wk_dir, genome_version))

	para_li = list(chr=chr_chars, QUAL_thresh=QUAL_threshs, maxIS_Mb=maxIS_Mbs)
	quality_tab = data.table(do.call(expand.grid, para_li))
	for(i in 1:ncol(quality_tab)) quality_tab[[i]] = as.vector(quality_tab[[i]])

	quality_tab$para = paste0(quality_tab$QUAL_thresh, '_', quality_tab$maxIS_Mb)

	####################################

	block_ratio_thresh = 0.01 ## it was 0.05 before 02.08.2021
	switch_err_thresh = 0.05 ## not sure whether this is too big; previously it was set as 0.01
	relative_phased_thresh = 0.7
	chr_char = paste0('chr', chr_num)
	quality_tab = quality_tab[chr==chr_char]
	n = nrow(quality_tab)

	#################################### compute consistency. Some parameters the overall selected block can still be low, does not matter

	para_li = list(i=1:n, j=1:n)
	para_tab = data.table(do.call(expand.grid, para_li))
	para_tab = para_tab[i>j]

	vcf_tab_li = foreach(i=1:n) %dopar%
	{
		QUAL_i = quality_tab$QUAL_thresh[i]
		max_IS_i = quality_tab$maxIS_Mb[i]
		vcf_file = sprintf('%s/mega/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.hap.phased.VCF', wk_dir, chr_num, QUAL_i, max_IS_i)
		vcf_tab = try(vcf2tab(block_ratio_thresh=block_ratio_thresh, vcf_file, which_blocks=NA, gz=FALSE, phased=TRUE, save2wig=TRUE, consider_ps=TRUE, wig_ground_truth=NA, format_shapeit2_hap=FALSE))
		# if(class(vcf_tab) == "try-error")  vcf_tab = data.table(1) ## in case when all block ratios are small
		vcf_tab
	}

	names(vcf_tab_li) = paste0('p', quality_tab$para)

	####################################

	relative_phased = sapply(vcf_tab_li, nrow) / max(sapply(vcf_tab_li, nrow))
	n_block = sapply(vcf_tab_li, function(v) length(unique(v$PS)))
	quality_tab$relative_phased = 0
	quality_tab$relative_phased = relative_phased
	quality_tab$n_block = n_block

	####################################

	mean_switch_mat = mean_identical_mat = mean_identical_mat_switched = matrix(0, n, n)

	identical_vals = foreach(k=1:nrow(para_tab), .combine=rbind) %dopar%
	{
		i = para_tab[[1]][k]
		j = para_tab[[2]][k]
		vcf_tab_x = vcf_tab_li[[i]][, mget(c('POS', 'hap_val'))]
		vcf_tab_y = vcf_tab_li[[j]][, mget(c('POS', 'hap_val'))]
		vcf_tab_xy = merge(vcf_tab_x, vcf_tab_y, by='POS')[order(POS)]
		identical = mean(sign(vcf_tab_xy[[3]])==sign(vcf_tab_xy[[2]]))
		t(c(i, j, max(identical, 1 - identical)))
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
		t(c(i, j, mean(abs(diff(vcf_tab_diff))))) ## maybe to redefine this value using the official switch error
	}

	mean_switch_mat[switch_vals[,1:2]] = switch_vals[,3]
	mean_switch_mat = as.matrix(my_forceSymmetric(mean_switch_mat))

	rownames(mean_switch_mat) = colnames(mean_switch_mat) = rownames(mean_identical_mat) = colnames(mean_identical_mat) = paste0(quality_tab$QUAL_thresh, '_', quality_tab$maxIS_Mb)

	dist_fun = function(A) return(as.dist(A))
	# branches = cutree(hclust(dist_fun(mean_switch_mat), method='average'), k=2)
	# branch_1 = names(which(branches==1))
	# branch_2 = names(which(branches==2))
	para_max = quality_tab[relative_phased==1]$para[1] ## the largest amount of SNP phased
	# branch2choose = ifelse( para_max %in% branch_1, 1, 2 )
	# # branch2choose = which.min(c(mean(mean_switch_mat[branch_1, branch_1]), mean(mean_switch_mat[branch_2, branch_2])))
	# para2choose = names( which(branches == branch2choose) )

	candidate_low_switch = names(which(mean_switch_mat[para_max, ] < switch_err_thresh))
	quality_tab$flag = "candidate"
	quality_tab[!(para %in% candidate_low_switch)]$flag = "err"

	# ####################################

	# quality_tab$switch_err = unname(apply(mean_switch_mat, 1, function(v) mean(head(sort(v[v>0]), 1)))) ## should be at least one is consistent with it
	# quality_tab$flag = 'candidate'
	# quality_tab[switch_err > switch_err_thresh]$flag = 'err'
	# # quality_tab[!(para %in% para2choose) ]$flag = 'err'
	# quality_tab[relative_phased < 0.9]$flag = 'err'

	# if(all(quality_tab$flag=='err')) quality_tab[order(switch_err)][1:2, 'flag'] = 'candidate' ## choose at least 2 candidate

	quality_tab[relative_phased < relative_phased_thresh]$flag = paste0(quality_tab[relative_phased < relative_phased_thresh]$flag, '&short')
	# if(all(quality_tab$flag=='err')) quality_tab[order(-relative_phased)][1:2, 'flag'] = 'candidate' ## choose at least 2 candidate

	if(sum(quality_tab$flag=='candidate')==1) quality_tab[order(-relative_phased)][1:2, 'flag'] = 'candidate' ## choose at least 2 candidate

	# keep_indices = which(quality_tab$n_block > 1)

	cellnote = matrix(substr(mean_switch_mat*100, 1, 3), nrow(mean_switch_mat))
	pdf(consist_plot_f, width=9, height=9)
	x = (mean_switch_mat^(1/3))
	rownames(x) = paste0(rownames(x), '( ', gsub('err', 'err', quality_tab$flag))
	try(heatmap.2(notecex=0.6, margins=c(12,12), revC=TRUE, key=0, hclustfun=my_hclust, distfun=dist_fun, cellnote=cellnote, x=x, dendrogram="both", col="bluered", trace="none", main="", ylab="", xlab=""))
	dev.off()

	saveRDS(mean_switch_mat, file=consist_rds_f)

	#################################### ## err is excluded for computing fragile site

	fragile_info = get_fragile_sites(vcf_tab_li[which(quality_tab$flag=='candidate')], save_file=fragile_plot_f)
	fragile_sites = c(0, fragile_info$fragile_sites, Inf)

	## when indeed n_fragile is 0
	if(any(is.na(fragile_sites)) | length(fragile_sites)==2) fragile_sites = c(0, 50, Inf)
	fragile_sites = sort(fragile_sites)

	n_fragile = length(fragile_sites) - 2
	vcf_tab_max = vcf_tab_li[[paste0('p', para_max)]]

	#####################################

	quality_tab = quality_tab[order(-relative_phased)]
	flip_info_li = list() 
	# for(index in 1:1) ## only choose vcf_tab_max to do flip
	{
	  vcf_tab = vcf_tab_max
	  vcf_tab$hap_val_bp = vcf_tab$hap_val

	  vcf_tab_bp = vcf_tab
	  keep_pos = vcf_tab$POS
	  QUAL = quality_tab$QUAL_thresh[1]
	  max_IS = quality_tab$maxIS_Mb[1]

	  switch = 0
	  switch_pos = NULL
	  switch_vec = NULL
	  stat_12 = list()

	  for(f in 1:n_fragile)
	  {

	  	  keep_range = c(fragile_sites[f], fragile_sites[f+2])
	  	  switch_range = c(fragile_sites[f+1], fragile_sites[f+2])

		  # snps_used_1 = generate_diploid_mat_pat_info(which_blocks=1:100, wk_dir, chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA, keep_pos=keep_pos, keep_range=keep_range)
		  snps_used_1 = generate_diploid_mat_pat_info_phased_vcf_tab(vcf_tab, file.path(wk_dir, 'mega'), chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA, keep_pos=keep_pos, keep_range=keep_range)
		  command_diploid = sprintf('%s %s %s %s %s %s', my_diploid_helper, wk_dir, chr_num, QUAL, max_IS, repo_dir)

		  system(command_diploid)
		  stat_1 = return_stat(wk_dir, chr_num, QUAL, max_IS)

		  # snps_used_2 = generate_diploid_mat_pat_info(which_blocks=1:100, wk_dir, chr_num, QUAL, max_IS, switch_range=switch_range, switch_pos=NA, keep_pos=keep_pos, keep_range=keep_range)
		  snps_used_2 = generate_diploid_mat_pat_info_phased_vcf_tab(vcf_tab, file.path(wk_dir, 'mega'), chr_num, QUAL, max_IS, switch_range=switch_range, switch_pos=NA, keep_pos=keep_pos, keep_range=keep_range)
		  command_diploid = sprintf('%s %s %s %s %s %s', my_diploid_helper, wk_dir, chr_num, QUAL, max_IS, repo_dir)
		  system(command_diploid)
		  stat_2 = return_stat(wk_dir, chr_num, QUAL, max_IS)

		  stat_12[[f]] = rbind(stat_1, stat_2)
		  switch = 1*(stat_1$Median > stat_2$Median)
		  switch_vec = c(switch_vec, switch)
	  }
	  
	  switch_vec_final = switch_vec
	  if(n_fragile > 1) {for(z in 1:(n_fragile-1)) if(switch_vec_final[z]==1) switch_vec_final[z+1] = abs(1-switch_vec_final[z+1])}
	  
	  for(z in 1:n_fragile) ## vcf_tab has been switched here
	  {
	  	switch_range = c(fragile_sites[z+1], fragile_sites[z+2])*1E6
	  	if(switch_vec_final[z]==1) vcf_tab[between(POS, switch_range[1], switch_range[2])]$hap_val = -1*vcf_tab[between(POS, switch_range[1], switch_range[2])]$hap_val
	  }

	  ## compute stat for the switched vcf_stat
	  snps_switch = generate_diploid_mat_pat_info_phased_vcf_tab(is_opt=TRUE, vcf_tab, file.path(wk_dir, 'mega'), chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA, keep_pos=keep_pos)
	  command_diploid = sprintf('%s %s %s %s %s %s', my_diploid_helper, wk_dir, chr_num, QUAL, max_IS, repo_dir)
	  system(command_diploid)
	  stat_switch = return_stat(wk_dir, chr_num, QUAL, max_IS)

	  ## compute stat for the original vcf_stat

	  stat_ori = stat_switch
	  if(any(switch_vec_final!=0))
	  {
		  snps_used_ori = generate_diploid_mat_pat_info_phased_vcf_tab(vcf_tab_bp, file.path(wk_dir, 'mega'), chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA, keep_pos=keep_pos)
		  command_diploid = sprintf('%s %s %s %s %s %s', my_diploid_helper, wk_dir, chr_num, QUAL, max_IS, repo_dir)
		  system(command_diploid)
		  stat_ori = return_stat(wk_dir, chr_num, QUAL, max_IS)
	  }

	  info = list(stat12=stat_12, stat_switch=stat_switch, stat_ori=stat_ori, switch_vec=switch_vec, vcf_tab=vcf_tab)
	  flip_info_li[[1]] = info
	}

	names(flip_info_li) = para_max
	flip_info_li[['fragile_sites']] = fragile_sites
	save(flip_info_li, file=fragile_flip_f)

