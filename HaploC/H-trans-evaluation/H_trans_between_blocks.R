    #!/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    wk_dir = args[1]
    ID = basename(wk_dir)

    ## h-trans profile between haplo-blocks

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
    source(sprintf('%s/header_v2.R', HaploC_dir))

	####################################

	bin_size = 100E3
	identity = c('QUAL=10_maxIS=50', 'opt')[1]

	wk_dir_mega = file.path(wk_dir, 'mega')
    htrans_stat_file = sprintf('%s/plots/htrans_ratio/htrans_ratio_%s.Rdata', wk_dir_mega, identity)

    ####################################

    load(htrans_stat_file)

    ####################################

    fragile_flip_f = sprintf('%s/mega/HapCut/fragile_sites/fragile_sites_flip_chr%s.Rdata', wk_dir, chr_num)
    load(fragile_flip_f)

    break_points = unname(flip_info_li$fragile_sites)
    break_bins = floor(break_points[between(break_points, 1, 1E7)]*1E6 / bin_size)
    break_bins = c(break_bins, 1E7)

    binned_chr = htrans_stat$chr1.binned
    tab2plot = binned_chr$diff
    # tab2plot = tab2plot[bin_index_x!=bin_index_y]

    tab2plot$block_x <- findInterval(tab2plot$bin_index_x, break_bins) + 1
    tab2plot$block_y <- findInterval(tab2plot$bin_index_y, break_bins) + 1
	tab2plot$block_xy = apply(tab2plot[, c('block_x', 'block_y')], 1, function(v) paste0(sort(v), collapse=''))

    # tab2plot = tab2plot[block_x!=block_y]
    tab2plot[, n_block_xy:=sum(n_contact), by='block_xy']
    tab2plot = unique(tab2plot, by='block_xy')

    n = length(break_points) - 1
    mat = matrix(0, n, n)
    for(i in 1:n)
    for(j in 1:n)
    {	
    	val = tab2plot[block_x==i][block_y==j]$n_block_xy
    	if(length(val)==0) next
    	mat[i,j] = mat[j,i] = val
    }

    rownames(mat) = colnames(mat) = round(break_points[2:(n+1)], 2)

    htrans_plot_file = sprintf('%s/plots/htrans_ratio/htrans_count_%s.pdf', wk_dir_mega, identity)
    pdf(htrans_plot_file, width=9, height=9)
   	try(heatmap.2(Rowv=FALSE, Colv=FALSE, notecex=1, margins=c(12,12), revC=TRUE, key=0, cellnote=mat, x=mat^(1/5), dendrogram="none", col="bluered", trace="none", main=identity, ylab="", xlab=""))
   	dev.off()

