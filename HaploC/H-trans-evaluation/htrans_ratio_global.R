    # !/usr/bin/env Rscript
    args = commandArgs(trailingOnly=TRUE)
    ## This function plots the h_trans ratio as a function of distance

    ###################################

    wk_dir = args[1]
    identity = args[2]
    chr_num = args[3]
    which_fun = args[4]

    chr_name = paste0("chr", chr_num)
    chr_chars = paste0('chr', c(1:22, 'X'))

    ###################################

    try(silent=TRUE, library(doParallel, quietly = T))
    try(silent=TRUE, library(data.table, quietly = T))
    registerDoParallel(cores=3)

    ################################

    identities = c("maternal_both", "paternal_both", "diff") ## both ends are maternal or patternal, or mat x pat
    bin_size = 1000E3
    mb_thresh = 80E6 ## aggregate values above 80Mb

    ###################################

    wk_dir_mega = file.path(wk_dir, 'mega')
    dir.create(file.path(wk_dir_mega, "plots/htrans_ratio/meta_data"), recursive=TRUE, showWarnings=FALSE)

    ###################################

    if(which_fun=="generate_data")
    {
        htrans_stat_file = sprintf('%s/plots/htrans_ratio/meta_data/htrans_ratio_%s.%s.rds', wk_dir_mega, identity, chr_name)
        if(!file.exists(htrans_stat_file))
        {
            htrans_stat = list()
            saveRDS(htrans_stat, file=htrans_stat_file)
        }

        freq_li = foreach(id=identities) %dopar%
        {
            contacts_tab = fread(sprintf("%s/HapCut/hic_phased/hic_chr%s_%s.%s.txt", wk_dir_mega, chr_num, identity, id))[, mget(c("V3", "V7"))]
            colnames(contacts_tab) = c('x', 'y')
            contacts_tab$dist = ceiling(abs(contacts_tab[[1]] - contacts_tab[[2]]) / bin_size)
            contacts_tab[dist > (mb_thresh/bin_size)]$dist = (mb_thresh/bin_size)
            freq = table(contacts_tab$dist)
            freq
        }
        names(freq_li) = identities

        ranges = range(as.numeric(Reduce(union, sapply(freq_li, names))))
        bins = as.character(ranges[1]:ranges[2])

        freq_tab = foreach(id=identities, .combine=cbind) %dopar%
        {
            unname(freq_li[[id]][bins])
        }

        freq_tab = data.table(freq_tab)
        colnames(freq_tab) = identities
        freq_tab[is.na(freq_tab)] = 0

        cis_trans_tab = data.table(cis=freq_tab$maternal_both + freq_tab$paternal_both, trans=freq_tab$diff, dist=bins)
        cis_trans_tab = cis_trans_tab[dist >0]

        cis_trans_tab$cis_cum = cumsum(cis_trans_tab$cis)
        cis_trans_tab$trans_cum = cumsum(cis_trans_tab$trans)

        cis_trans_tab$h_trans = cis_trans_tab$trans / (cis_trans_tab$cis + cis_trans_tab$trans)
        cis_trans_tab$h_trans_cum = cis_trans_tab$trans_cum / (cis_trans_tab$cis_cum + cis_trans_tab$trans_cum)

        htrans_stat = readRDS(htrans_stat_file)
        htrans_stat[[chr_name]] = cis_trans_tab
        saveRDS(htrans_stat, file=htrans_stat_file)

        ################################### between fragile sites and bins (100kb)

        fragile_flip_f = sprintf('%s/mega/HapCut/fragile_sites/fragile_sites_flip_chr%s.Rdata', wk_dir, chr_num)
        load(fragile_flip_f)

        break_points = sort(unname(flip_info_li$fragile_sites)*1E6)

        block_freq_li = foreach(id=identities) %dopar%
        {
            contacts_tab = fread(sprintf("%s/HapCut/hic_phased/hic_chr%s_%s.%s.txt", wk_dir_mega, chr_num, identity, id))[, mget(c("V3", "V7"))]
            colnames(contacts_tab) = c('x', 'y')
            contacts_tab[y <= x, c("x", "y") := .(y, x)]

            contacts_tab$block_x <- findInterval(contacts_tab$x, break_points)
            contacts_tab$block_y <- findInterval(contacts_tab$y, break_points)
            contacts_tab$bin_index_x = ceiling(contacts_tab$x / 100E3)
            contacts_tab$bin_index_y = ceiling(contacts_tab$y / 100E3)
            contacts_tab[, n_contact_bbbb:=length(x), by=c('block_x', 'block_y', 'bin_index_x', 'bin_index_y')] ## bbbb: block + bin | x and y
            bbbb_tab = unique(contacts_tab, by=c('block_x', 'block_y', 'bin_index_x', 'bin_index_y'))[, -c(1,2)]
            bbbb_tab[, n_contact_bb:=sum(n_contact_bbbb), by=c('block_x', 'block_y')] ## bb block_x and block_y
            bb_tab = unique(bbbb_tab, by=c('block_x', 'block_y'))[, c('block_x', 'block_y', 'n_contact_bb')]
            li = list(bbbb_tab=bbbb_tab, bb_tab=bb_tab)
            li
        }

        names(block_freq_li) = identities
        htrans_stat[[sprintf('%s.binned', chr_name)]] = block_freq_li
        saveRDS(htrans_stat, file=htrans_stat_file)
    }

    ################################

    if(which_fun=="generate_plot")
    {
        htrans_stat = foreach(chr_name=chr_chars, .combine=c) %do%
        {
            htrans_stat_file = sprintf('%s/plots/htrans_ratio/meta_data/htrans_ratio_%s.%s.rds', wk_dir_mega, identity, chr_name)
            if(!file.exists(htrans_stat_file)) return(NULL)
            htrans_stat = try(readRDS(htrans_stat_file))
            if(class(htrans_stat)[1]=='try-error') return(NULL)
            htrans_stat
        }

        htrans_plot_file = sprintf('%s/plots/htrans_ratio/htrans_ratio_%s.pdf', wk_dir_mega, identity)
        
        pdf(htrans_plot_file, width=10, height=4)
        par(mfrow=c(1, 3))

        for(chr_name in chr_chars)
        {
            if(!(chr_name %in% names(htrans_stat))) next

            cis_trans_tab = htrans_stat[[chr_name]]
            xlab = sprintf("%s | bin_size:%sMb", chr_name, bin_size/1E6)
            plot(main="H-trans at a given Mb distance", cis_trans_tab$dist, cis_trans_tab$h_trans, type="l", col="red", xlab=xlab, ylab="htrans_ratio", ylim=c(0, 0.5), log='x')
            abline(v=30, h=0.05, col="blue", lty=2)
            legend("topleft", legend=chr_name, bty="n")

            plot(main="cumulative=T", cis_trans_tab$dist, cis_trans_tab$h_trans_cum, type="l", col="red", xlab=xlab, ylab="htrans_ratio", ylim=c(0, 0.2), log='x')
            abline(v=30, h=0.05, col="blue", lty=2)    
            legend("topleft", legend=chr_name, bty="n")

            htrans_ratio_global = tail(cis_trans_tab$trans_cum / (cis_trans_tab$trans_cum + cis_trans_tab$cis_cum), 1)
            barplot(htrans_ratio_global, ylim=c(0,0.3), yaxs="i", xaxs="i", col="red", space=c(1,0.5), main="global_htrans_ratio", xlim=c(1,1.5), xlab='global H-trans')
            legend("topleft", legend=format(htrans_ratio_global, digits=4), bty="n")
        }

        dev.off()
    }

    