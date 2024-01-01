
    ########################################################

    source(sprintf('%s/header/shared_header.R', repo_dir))
    options(scipen=20) ## this is required for the maxIS parameter for HaploC 

    ########################################################

    bamtools = "bamtools"
    samtools = "samtools"
    bcftools = 'bcftools'
    seqkit =  file.path(HaploC_dir, "/dependancies/seqkit/seqkit") ## download from here: https://bioinf.shenwei.me/seqkit/download/
    juicer_jar = file.path(juicerDir, "scripts/juicer_tools_1.19.02.jar")  

    HapCUT22  = 'hapcut2' ## conda version
    EXTRACTHAIRS = "extractHAIRS"  

    ########################################################

    free_exe = file.path(HaploC_dir, "/dependancies/freebayes") ## v1.3.6
    SHAPEIT2     = file.path(HaploC_dir, '/HapCUT2/SHAPEIT2.sh')
    HiC_repair_py = file.path(HaploC_dir, '/HapCUT2/my_hic_repair.py')
    my_diploid_helper = file.path(HaploC_dir, '/create_phased_maps/my_diploid_helper.sh')
    generate_mated_contacts = file.path(HaploC_dir, '/H-trans-evaluation/generate_mated_contacts.sh')
    intersect_merged_nodup = file.path(HaploC_dir, '/SNP/intersect_with_snp.sh')
    htrans_cis_info_global_fun = file.path(HaploC_dir, "/H-trans-evaluation/htrans_ratio_global.R")
    diploid_opt_intra_sh =  file.path(HaploC_dir, "/create_phased_maps/diploid_opt_intra_helper.sh")
    generate_random_mixed_sh = file.path(HaploC_dir, "/create_phased_maps/generate_random_mixed_hic.sh")

    #########################################################

    get_fragile_sites = function(vcf_tab_li, save_file)
    {
        my_merge = function(...) merge(..., by='POS', all=TRUE)
        vcf_tab_li_slim = vcf_tab_li
        len = length(vcf_tab_li_slim)
        n_col = len + 1

        for(i in 1:len)
        {
            n_blocks = length(unique(vcf_tab_li[[i]]$block_index))

            set.seed(i) ## for different vcfs
            rand_vals = abs(rnorm(n_blocks, 10, 10)) ## add a random value for the same block
            vcf_tab_li_slim[[i]]$hap_val = (vcf_tab_li_slim[[i]]$hap_val)*(rand_vals[(vcf_tab_li_slim[[i]]$block_index)]) ## sam abs value if soming from the same hap block
            vcf_tab_li_slim[[i]] = vcf_tab_li_slim[[i]][, mget(c('POS', 'hap_val'))]
            colnames(vcf_tab_li_slim[[i]])[2] = paste0(colnames(vcf_tab_li_slim[[i]])[2], '_', names(vcf_tab_li)[i])
        }

        vcf_combined = Reduce(my_merge, vcf_tab_li_slim)
        vcf_combined = vcf_combined[order(POS)]

        # signs = sign(vcf_combined$hap_val_p1_1000)

        na_count = apply(vcf_combined[,-1], 2, function(v) sum(is.na(v)) )
        low_na_para = names(which.min(na_count))

        signs = unname(unlist(sign(vcf_combined[, mget(low_na_para)])))
        for(i in 2:ncol(vcf_combined)) vcf_combined[[i]] = signs*vcf_combined[[i]]


        vcf_combined$ALL =  apply(abs(vcf_combined[, -c(1, 2), with=FALSE] - vcf_combined[, -c(1, n_col), with = FALSE]), 1, mean, na.rm=TRUE)
        vcf_combined = vcf_combined[!is.nan(ALL)]
        vcf_combined[ALL<1]$ALL = vcf_combined[ALL<1]$ALL + 1


        vcf_combined[, `:=` (segment_ratio = .N/nrow(vcf_combined)), by = ALL]
        vcf_combined = vcf_combined[segment_ratio > 0.05]
        vcf_combined$segment_ratio = substr(vcf_combined$segment_ratio, 1, 10)

        centoids_raw = tapply(vcf_combined$POS/1E6, vcf_combined$segment_ratio, median)
        n_segments = length(centoids_raw)
        segments = 1:n_segments
        names(segments) = names(sort(centoids_raw))
        vcf_combined$segment_index = unname(segments[vcf_combined$segment_ratio])
        centoids = tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, median)
        fivenums = do.call(rbind, tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, summary))

        ten_percentile = c(sort(tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, quantile, 0.01)*1E6), Inf)
        uniq_ALLs = unique(vcf_combined$ALL) ## already ordered based on first appearance of POS


        vcf_combined$stair_val = 0
        for(z in 1:length(uniq_ALLs)) vcf_combined[between(POS, ten_percentile[z], ten_percentile[z+1])]$stair_val = apply(vcf_combined[between(POS, ten_percentile[z], ten_percentile[z+1]), 'ALL'], 1, function(v) sum(setdiff(uniq_ALLs[1:z], v)))
        
        # vcf_combined$stair_val = substr(vcf_combined$stair_val, 1, 20)
        vcf_combined[, `:=` (new_count = .N/nrow(vcf_combined)), by = stair_val]
        vcf_combined = vcf_combined[new_count > 0.05] ## ignore those with less than 5%
        vcf_combined$segment_index = match(vcf_combined$stair_val, unique(vcf_combined$stair_val))


        n_segments = length(unique(vcf_combined$segment_index))
        fragile_sites = NULL
        for(s in 1:(n_segments-1)) 
        {
            quant_1 = quantile(vcf_combined[segment_index==s]$POS, 0.95) 
            quant_2 = quantile(vcf_combined[segment_index==(s+1)]$POS, 0.05)
            mid = (quant_1 + quant_2)/2/1E6
            fragile_sites = c(fragile_sites, mid)
        }

        fragile_info = list(fragile_sites=fragile_sites, vcf_combined=vcf_combined)

        pdf(save_file, width=10, height=6)
        plot(xaxt="n", xlab='Mb', ylab='segment pseudo val', vcf_combined$POS/1E6, vcf_combined$ALL, pch='.', col='red', cex=0.1, type='p')
        axis(1, at=(1:100)*5, las=2)
        abline(v=fragile_sites, lty=2, col='blue')
        # abline(v=(1:100)*5, lty=2, col=grey(0.9))
        # abline(v=ten_percentile/1E6, lty=2, col='orange')
        # lines(vcf_combined$POS/1E6, vcf_combined$stair_val, pch='.', col='blue', cex=0.1, type='p')
        graphics.off()

        return(fragile_info)
    }

    get_fragile_sites_simple_case = function(vcf_tab_li, save_file)
    {
        my_merge = function(...) merge(..., by='POS', all=TRUE)
        vcf_tab_li_slim = vcf_tab_li
        len = length(vcf_tab_li_slim)
        n_col = len + 1

        for(i in 1:len)
        {
            n_blocks = length(unique(vcf_tab_li[[i]]$block_index))

            set.seed(i) ## for different vcfs
            rand_vals = abs(rnorm(n_blocks, 10, 10)) ## add a random value for the same block
            vcf_tab_li_slim[[i]]$hap_val = (vcf_tab_li_slim[[i]]$hap_val)*(rand_vals[(vcf_tab_li_slim[[i]]$block_index)])
            vcf_tab_li_slim[[i]] = vcf_tab_li_slim[[i]][, mget(c('POS', 'hap_val'))]
            colnames(vcf_tab_li_slim[[i]])[2] = paste0(colnames(vcf_tab_li_slim[[i]])[2], '_', names(vcf_tab_li)[i])
        }

        vcf_combined = Reduce(my_merge, vcf_tab_li_slim)
        vcf_combined = vcf_combined[order(POS)]

        vcf_combined$ALL =  apply(abs(vcf_combined[, -c(1, 2), with=FALSE] - vcf_combined[, -c(1, n_col), with = FALSE]), 1, mean, na.rm=TRUE)
        vcf_combined = vcf_combined[!is.nan(ALL)]
        vcf_combined[ALL<1]$ALL = vcf_combined[ALL<1]$ALL + 1


        vcf_combined[, `:=` (segment_ratio = .N/nrow(vcf_combined)), by = ALL]
        # vcf_combined = vcf_combined[segment_ratio > 0.05]
        # segs = segment(vcf_combined$ALL, maxseg=6, maxk=15)


        vcf_combined$segment_ratio = substr(vcf_combined$segment_ratio, 1, 10)

        centoids_raw = tapply(vcf_combined$POS/1E6, vcf_combined$segment_ratio, median)
        n_segments = length(centoids_raw)
        segments = 1:n_segments
        names(segments) = names(sort(centoids_raw))
        vcf_combined$segment_index = unname(segments[vcf_combined$segment_ratio])
        centoids = tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, median)
        # fivenums = do.call(rbind, tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, summary))

        # ten_percentile = c(sort(tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, quantile, 0.1)*1E6), Inf)
        # uniq_ALLs = unique(vcf_combined$ALL) ## already ordered based on first appearance of POS


        # vcf_combined$stair_val = 0
        # for(z in 1:length(uniq_ALLs)) vcf_combined[between(POS, ten_percentile[z], ten_percentile[z+1])]$stair_val = apply(vcf_combined[between(POS, ten_percentile[z], ten_percentile[z+1]), 'ALL'], 1, function(v) sum(setdiff(uniq_ALLs[1:z], v)))
        # vcf_combined$stair_val = substr(vcf_combined$stair_val, 1, 3)
        # vcf_combined$segment_index = match(vcf_combined$stair_val, unique(vcf_combined$stair_val))

        # n_segments = length(unique(vcf_combined$segment_index))
        fragile_sites = NULL
        for(s in 1:(n_segments-1)) 
        {
            quant_1 = quantile(vcf_combined[segment_index==s]$POS, 0.95) 
            quant_2 = quantile(vcf_combined[segment_index==(s+1)]$POS, 0.05)
            mid = (quant_1 + quant_2)/2
            fragile_sites = c(fragile_sites, mid)
        }

        fragile_info = list(fragile_sites=fragile_sites, vcf_combined=vcf_combined)

        pdf(save_file, width=10, height=6)
        plot(xaxt="n", xlab='Mb', ylab='segment pseudo val', vcf_combined$POS/1E6, vcf_combined$ALL, pch='.', col='red', cex=0.1, type='p')
        axis(1, at=(1:100)*5, las=2)
        abline(v=fragile_sites/1E6, lty=2, col='blue')
        # abline(v=(1:100)*5, lty=2, col=grey(0.9))
        graphics.off()

        return(fragile_info)
    }


    ma <- function(x, n = 5) {filter(x, rep(1 / n, n), sides = 2)} ## https://stackoverflow.com/questions/743812/calculating-moving-average


    return_stat = function(wk_dir, chr_num, QUAL, max_IS)
    {
      stat_file = sprintf('%s/mega/HapCut/hic_phased/hic_chr%s_QUAL=%s_maxIS=%s.stat.txt', wk_dir, chr_num, QUAL, max_IS)
      stats_tab = data.table(t((as.numeric(t(fread(stat_file))[2,]))))
      colnames(stats_tab) = c('mat', 'pat', 'diff')
      stats_tab$mat_pat = stats_tab$mat + stats_tab$pat
      stats_tab$diff_ratio = stats_tab$diff / (stats_tab$mat_pat + stats_tab$diff)
      stats_tab

      diff_file = sprintf('%s/mega/HapCut/hic_phased/hic_chr%s_QUAL=%s_maxIS=%s.diff.txt', wk_dir, chr_num, QUAL, max_IS)
      locus_dist = fread(diff_file)
      summary_stat = summary(abs(locus_dist[[3]] - locus_dist[[7]])) / 1E6
      cbind(stats_tab, setDT(as.list(summary_stat[2:4]))[])
    }

    generate_diploid_mat_pat_info = function(wk_dir, chr_num, switch_range=NA, switch_pos=NA, QUAL, max_IS, keep_pos=NA, keep_range=NA, which_blocks=1, add_nphased=FALSE) ## if I want to switch the haplotype position within a given range
    {
        switch_range = switch_range*1E6 ## Mb
        keep_range = keep_range*1E6
        chr_num = chr_num = sort(chr_num) ## for a single chr
        hap_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.hap', wk_dir, chr_num, QUAL, max_IS)
        vcf_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.hap.phased.VCF', wk_dir, chr_num, QUAL, max_IS)
        vcf_tab = vcf2tab(vcf_file, QUAL_thresh=10, gz=FALSE, phased=FALSE, save2wig=FALSE, consider_ps=FALSE, wig_ground_truth=NA, format_shapeit2_hap=FALSE)

        SNP_chrs = foreach(chr_num=chr_num, .combine=rbind) %dopar%
        {
          # hap_file = sprintf('%s/phased_hap_intg/hic_chr%s.hap', wk_dir, chr_num)

          blocks = data.table(system(intern=TRUE, sprintf('cut -f1-7 %s', hap_file)))
          blocks = blocks[!grepl('[****]', V1)]
          blocks$is_start = apply(blocks, 1, function(v) grepl('BLOCK', v)*1)
          blocks$block_index = cumsum(blocks$is_start)
          block_li_raw = split(blocks, by='block_index')

          block_tab_li = foreach(i=1:length(block_li_raw)) %do%
          {
              block = block_li_raw[[i]]
              block_tab = data.table(t(apply(block[-1], 1, function(u) strsplit(u, '\t')[[1]])))[, -1]
                
              # block_tab_max = block_tab[block_index==1] 
              block_tab$HET_block = paste0(block_tab[[1]], '|', block_tab[[2]])
              block_tab = block_tab[HET_block %in% c('0|1', '1|0'), -c(1:2)]
              colnames(block_tab)[1:2] = c('chr', 'POS')
              block_tab$POS = as.numeric(block_tab$POS)
              block_tab
          }

          block_tab_li = block_tab_li[order(-sapply(block_tab_li, nrow))]
          block_tab_li = block_tab_li[sapply(block_tab_li, nrow)>0]

          for(index in 1:length(block_tab_li)) block_tab_li[[index]]$block_index = index ## ranked by block size
          first_n_block_size = substr(head(sapply(block_tab_li, nrow), 5), 1, 7)
          cat('First 5 block size: ', first_n_block_size, '\n')
          which_blocks = intersect(which_blocks, 1:length(block_tab_li))
          block_tab_max = do.call(rbind, block_tab_li[which_blocks])
          block_tab_max
        }

        if(length(keep_pos) > 1) SNP_chrs = SNP_chrs[POS %in% keep_pos]
        if(length(keep_range) > 1) SNP_chrs = SNP_chrs[between(POS, min(keep_range), max(keep_range))]

        n0 = nrow(SNP_chrs)
        # #########################################################

        chr_pos_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.chr_pos', wk_dir, chr_num, QUAL, max_IS)
        snp_pat_mat_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.pat_mat', wk_dir, chr_num, QUAL, max_IS)

        if(add_nphased)
        {
            chr_pos_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.chr_pos.add_nphased', wk_dir, chr_num, QUAL, max_IS)
            snp_pat_mat_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.pat_mat.add_nphased', wk_dir, chr_num, QUAL, max_IS)
            nphased_SNPs_raw = vcf_tab[!(POS %in% SNP_chrs$POS)][, mget(c('chr', 'POS', 'REF', 'ALT'))]
            nphased_SNPs = data.table(chr=nphased_SNPs_raw$chr, POS=nphased_SNPs_raw$POS, V6=nphased_SNPs_raw$REF, V7=nphased_SNPs_raw$ALT, HET_block='0|1', block_index=7777777)
            SNP_chrs = rbind(SNP_chrs, nphased_SNPs)
        }

        # cat(chr_pos_file, '\n')
        cat('#snps used for computing pat_mat (not add_nphased):', n0, '\n')
        cat('#snps used for computing pat_mat (add_nphased):', nrow(SNP_chrs), '\n')

        SNP_li = split(SNP_chrs, by='chr')
        for(i in 1:length(SNP_li))
        {
            SNP_ind = SNP_li[[i]][order(POS)]
            chr_pos = as.data.table(t(as.data.frame(c(SNP_ind$chr[1], SNP_ind$POS))))
            SNP_ind$mat = SNP_ind$pat = ''

            SNP_ind[HET_block=='0|1']$pat = SNP_ind[HET_block=='0|1']$V6
            SNP_ind[HET_block=='1|0']$pat = SNP_ind[HET_block=='1|0']$V7
            
            SNP_ind[HET_block=='0|1']$mat = SNP_ind[HET_block=='0|1']$V7
            SNP_ind[HET_block=='1|0']$mat = SNP_ind[HET_block=='1|0']$V6


            if(!is.na(switch_range[1]))
            {
                SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='0|1']$pat = SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='0|1']$V7
                SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='1|0']$pat = SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='1|0']$V6

                SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='0|1']$mat = SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='0|1']$V6
                SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='1|0']$mat = SNP_ind[between(POS, min(switch_range), max(switch_range))][HET_block=='1|0']$V7
            }

            if(is.numeric(switch_pos[1]))
            {
                cat('I switched some regions\n')
                SNP_ind[POS %in% switch_pos][HET_block=='0|1']$pat = SNP_ind[POS %in% switch_pos][HET_block=='0|1']$V7
                SNP_ind[POS %in% switch_pos][HET_block=='1|0']$pat = SNP_ind[POS %in% switch_pos][HET_block=='1|0']$V6

                SNP_ind[POS %in% switch_pos][HET_block=='0|1']$mat = SNP_ind[POS %in% switch_pos][HET_block=='0|1']$V6
                SNP_ind[POS %in% switch_pos][HET_block=='1|0']$mat = SNP_ind[POS %in% switch_pos][HET_block=='1|0']$V7
            }

            SNP_ind$chr_pos = paste0(SNP_ind$chr, ':', SNP_ind$POS)

            fwrite(SNP_ind[, mget(c('chr_pos', 'mat', 'pat'))], file=snp_pat_mat_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))        
            fwrite(chr_pos, file=chr_pos_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))
        }

        if(!add_nphased) return(SNP_ind[, mget(c('chr_pos', 'mat', 'pat'))])
        if(add_nphased) return(nphased_SNPs)
    }

    generate_diploid_mat_pat_info_n_phased = function(vcf_dir, chr_num) ## this will be used to extract merged_nodup.txt that overlaps snps
    {
        # vcf_file = sprintf('%s/VCFs/%s.fastq.gz.pos_sorted.chr%s.vcf.HET.bi.QUAL=10', wk_dir, prefix, chr_num)
        vcf_file = sprintf('%s/mega.chr%s.HET.bi.QUAL=10.vcf', vcf_dir, chr_num)
        vcf_tab = vcf2tab(vcf_file, QUAL_thresh=-1, gz=FALSE, phased=FALSE, save2wig=FALSE, consider_ps=FALSE, wig_ground_truth=NA, format_shapeit2_hap=FALSE)

        if(nrow(vcf_tab)==0) return(NULL)

        chr_pos_file = gsub('[.]vcf', '.chr_pos', vcf_file)
        snp_pat_mat_file = gsub('[.]vcf', '.pat_mat', vcf_file)

        # #########################################################
        
        vcf_tab$pat = 'N' ## Just pseudo value
        vcf_tab$mat = 'N' ## Just pseudo value
        vcf_tab$chr_pos = paste0(vcf_tab$chr, ':', vcf_tab$POS)

        keep_names = c('chr', 'POS', 'pat', 'mat', 'chr_pos')
        vcf_tab = vcf_tab[, mget(keep_names)]
        SNP_chrs = vcf_tab

        SNP_li = split(SNP_chrs, by='chr')
        for(i in 1:length(SNP_li))
        {
            SNP_ind = SNP_li[[i]][order(POS)]
            chr_pos = as.data.table(t(as.data.frame(c(SNP_ind$chr[1], SNP_ind$POS))))

            SNP_ind$chr_pos = paste0(SNP_ind$chr, ':', SNP_ind$POS)

            fwrite(SNP_ind[, mget(c('chr_pos', 'mat', 'pat'))], file=snp_pat_mat_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))        
            fwrite(chr_pos, file=chr_pos_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))
            # cat(i, '\n')
        }

        cat(sprintf('chr%s done\n', chr_num))
        return(c(chr_pos=chr_pos_file, snp_pat_mat=snp_pat_mat_file))
    }

    generate_diploid_mat_pat_info_phased_vcf_tab = function(phased_vcf_tab, wk_dir, chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA, keep_pos=NA, keep_range=NA, is_opt=FALSE, is_ground_truth=FALSE) ## if I want to switch the haplotype position within a given range
    {
        keep_range = keep_range*1E6
        switch_range = switch_range*1E6 ## Mb

        chr_pos_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.chr_pos', wk_dir, chr_num, QUAL, max_IS)
        snp_pat_mat_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.pat_mat', wk_dir, chr_num, QUAL, max_IS)

        if(is_opt)
        {
            chr_pos_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_opt.chr_pos', wk_dir, chr_num)
            snp_pat_mat_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_opt.pat_mat', wk_dir, chr_num)
        }

        if(is_ground_truth)
        {
            chr_pos_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_ground_truth.chr_pos', wk_dir, chr_num)
            snp_pat_mat_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_ground_truth.pat_mat', wk_dir, chr_num)
        }


        SNP_ind = phased_vcf_tab[!is.na(hap_val)]	

        SNP_ind$chr_pos = paste0(SNP_ind$chr, ':', SNP_ind$POS)
        if(!is.na(keep_pos[1])) SNP_ind = SNP_ind[POS %in% keep_pos]       
        if(!is.na(switch_range[1])) SNP_ind[between(POS, switch_range[1], switch_range[2])]$hap_val = -1*SNP_ind[between(POS, switch_range[1], switch_range[2])]$hap_val
        if(is.numeric(switch_pos[1])) SNP_ind[POS %in% switch_pos]$hap_val = -1*SNP_ind[POS %in% switch_pos]$hap_val
        if(length(keep_range) > 1) SNP_ind = SNP_ind[between(POS, min(keep_range), max(keep_range))]

        SNP_ind$mat = SNP_ind$pat = ''

        SNP_ind[hap_val>0]$pat = SNP_ind[hap_val>0]$REF
        SNP_ind[hap_val<0]$pat = SNP_ind[hap_val<0]$ALT
            
        SNP_ind[hap_val>0]$mat = SNP_ind[hap_val>0]$ALT
        SNP_ind[hap_val<0]$mat = SNP_ind[hap_val<0]$REF
        SNP_ind = SNP_ind[order(POS)]

        i = 1
        chr_pos = as.data.table(t(as.data.frame(c(SNP_ind$chr[1], SNP_ind$POS))))

        fwrite(SNP_ind[, mget(c('chr_pos', 'mat', 'pat'))], file=snp_pat_mat_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))        
        fwrite(chr_pos, file=chr_pos_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))
        
        return(SNP_ind[, mget(c('chr_pos', 'mat', 'pat'))])
    }

    generate_diploid_mat_pat_info_slim = function(phased_vcf_tab, wk_dir, chr_num, QUAL, max_IS, switch_range=NA, switch_pos=NA) ## if I want to switch the haplotype position within a given range
    {
        switch_range = switch_range*1E6 ## Mb
        chr_pos_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.chr_pos.switch', wk_dir, chr_num, QUAL, max_IS)
        snp_pat_mat_file = sprintf('%s/HapCut/phased_hap_intg/hic_chr%s_QUAL=%s_maxIS=%s.pat_mat.switch', wk_dir, chr_num, QUAL, max_IS)

        SNP_chrs = phased_vcf_tab
        SNP_li = split(SNP_chrs, by='chr')
        for(i in 1:length(SNP_li))
        {
            SNP_ind = SNP_li[[i]][order(POS)]
            chr_pos = as.data.table(t(as.data.frame(c(SNP_ind$chr[1], SNP_ind$POS))))
            SNP_ind$mat = SNP_ind$pat = ''

            SNP_ind[hap=='0|1']$pat = SNP_ind[hap=='0|1']$REF
            SNP_ind[hap=='1|0']$pat = SNP_ind[hap=='1|0']$ALT
            
            SNP_ind[hap=='0|1']$mat = SNP_ind[hap=='0|1']$ALT
            SNP_ind[hap=='1|0']$mat = SNP_ind[hap=='1|0']$REF


            if(!is.na(switch_range[1]))
            {
                SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='0|1']$pat = SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='0|1']$V7
                SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='1|0']$pat = SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='1|0']$V6

                SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='0|1']$mat = SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='0|1']$V6
                SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='1|0']$mat = SNP_ind[between(POS, min(switch_range), max(switch_range))][hap=='1|0']$V7
            }

            if(!is.na(switch_pos[1]))
            {
                SNP_ind[POS %in% switch_pos][hap=='0|1']$pat = SNP_ind[POS %in% switch_pos][hap=='0|1']$V7
                SNP_ind[POS %in% switch_pos][hap=='1|0']$pat = SNP_ind[POS %in% switch_pos][hap=='1|0']$V6

                SNP_ind[POS %in% switch_pos][hap=='0|1']$mat = SNP_ind[POS %in% switch_pos][hap=='0|1']$V6
                SNP_ind[POS %in% switch_pos][hap=='1|0']$mat = SNP_ind[POS %in% switch_pos][hap=='1|0']$V7
            }

            fwrite(SNP_ind[, mget(c('chr_pos', 'mat', 'pat'))], file=snp_pat_mat_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))        
            fwrite(chr_pos, file=chr_pos_file, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ', append=ifelse(i!=1, TRUE, FALSE))
        }
        return(SNP_chrs)
    }
    
    ## return vcf_tab or save the phased as wig
    ## block_ratio_thresh: blocks to select
    vcf2tab = function(vcf_file, gz=FALSE, phased=TRUE, save2wig=FALSE, consider_ps=FALSE, QUAL_thresh=-1, wig_ground_truth=NA, format_shapeit2_hap=FALSE, block_ratio_thresh=NA, which_blocks=1, vcf_tab_only=TRUE)
    {
        hap_val = c(1, -1, 0)
        names(hap_val) = c('0|1', '1|0', '0/1')

        if(format_shapeit2_hap==TRUE)
        {
            wig_file = gsub('.haps', '.wig', vcf_file, ignore.case=TRUE)
            vcf_tab = read.table(vcf_file, sep=' ', header=FALSE)
            vcf_tab = vcf_tab[, c(1,3:7)]
            colnames(vcf_tab)[c(1,2)] = c('chr', 'POS')
            vcf_tab$hap = paste0(vcf_tab[[5]], '|', vcf_tab[[6]]) 
            n_phased = nrow(vcf_tab)
            vcf_tab = vcf_tab[, c('chr', 'POS', 'hap')]
            vcf_tab = vcf_tab[which(vcf_tab$hap %in% c('0/1', '0|1', '1|0')), ] ## very small cases can be "./." 
            vcf_tab$hap_val = unname(hap_val[vcf_tab$hap])
            wig_tab = data.frame(chr=gsub('chrchr', 'chr', paste0('chr', vcf_tab[[1]])), pos_start=vcf_tab$POS - 1, pos_end=vcf_tab$POS, hap=vcf_tab$hap_val)
            wig_tab = na.omit(wig_tab)
            write.table(wig_tab, file=wig_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
            return(vcf_tab)
        }

        if(!format_shapeit2_hap)
        {
            if(gz==TRUE)
            {
                wig_file = gsub('.vcf.gz', '.wig', vcf_file, ignore.case=TRUE)
                my_skip = as.numeric(strsplit(system(intern=TRUE, sprintf('zgrep -n "#CHROM\tPOS" %s', vcf_file)), ':')[[1]][1]) - 1
            }

            if(gz==FALSE)
            {
                wig_file = gsub('[.]vcf', '.wig', vcf_file, ignore.case=TRUE)
                my_skip = as.numeric(strsplit(system(intern=TRUE, sprintf('grep -n "#CHROM\tPOS" %s', vcf_file)), ':')[[1]][1]) - 1
            }

            wig_diff_file = gsub('.vcf', '.diff.wig', vcf_file, ignore.case=TRUE)

            vcf_tab = data.table::fread(vcf_file, sep='\t', skip=my_skip)
            sample_name = tail(colnames(vcf_tab), 1)

            if(consider_ps)
            {
                PS = system(sprintf('bcftools query -f "[%%PS]\n" %s', vcf_file), intern=TRUE)
                vcf_tab$PS = PS
            }

            vcf_tab = vcf_tab[, mget(setdiff(colnames(vcf_tab), 'INFO'))]
            vcf_tab = vcf_tab[QUAL > QUAL_thresh]
        }

        vcf_tab$hap = substr(vcf_tab[[sample_name]], 1, 3)
        vcf_tab = vcf_tab[hap %in% c('0/1', '0|1', '1|0')] ## very small cases can be "./." 

        ## If wig_ground_truth is not NA
        if(length(wig_ground_truth) > 1) vcf_tab = vcf_tab[POS %in% wig_ground_truth$pos_end]

        colnames(vcf_tab)[c(1,2)] = c('chr', 'POS') 

        vcf_tab_input = vcf_tab
        n_total_input_SNPs =  nrow(vcf_tab_input) ## total number of variants in the VCF file input for phasing

        #########################################################

        if(consider_ps) 
        {
            # vcf_tab$PS = unname(sapply(vcf_tab[[sample_name]], function(v) tail(strsplit(v, ':')[[1]],1)))
            # vcf_tab[hap=='0/1']$PS = '.'
            vcf_tab = vcf_tab[PS!='.']

            vcf_tab[, `:=` (block_ratio = .N/nrow(vcf_tab)), by = PS]

            # vcf_tab = vcf_tab[order(-block_ratio, POS)]
            # vcf_tab$block_index = c(1, cumsum(diff(vcf_tab$block_ratio)!=0) + 1) ## this is a bug until 01-07-2021: because block with the same size with get the same block index, not good...
            vcf_tab = vcf_tab[order(PS, -block_ratio)]
            vcf_tab$block_index = c(1, cumsum(diff(as.numeric(vcf_tab$PS))!=0) + 1)
            ## should sort by block size

            block_sizes_order = 1:length(unique(vcf_tab$PS))
            names(block_sizes_order) = as.numeric(names(sort(table(vcf_tab$PS), decreasing=TRUE)))
            vcf_tab$block_index = unname(block_sizes_order[as.character(vcf_tab$PS)])

            ## there can also be blocks with only one SNP, because of the following situation:
          
            #             BLOCK: offset: 2471 len: 1154 phased: 3 SPAN: 12904954 fragments 2
            # 2471    -   -   chrX    24359592    A   C   0/1:40:33,7:33:1134:7:267:-12.3493,0,-90.2698   0   .   6.38    1
            # 2554    0   1   chrX    25011374    C   T   0/1:36:22,14:22:808:14:528:-37.0195,0,-62.1856  0   .   10.37   2
            # 3624    -   -   chrX    37264546    T   C   0/1:37:31,6:31:1187:6:234:-10.291,0,-95.9775    0   .   5.97    1
            # ******** 

            vcf_tab[, `:=` (block_range = paste0(substr(min(POS)/1E6, 1, 4), '_', substr(max(POS)/1E6, 1, 4))), by = PS]

            n_phased = nrow(vcf_tab)
            # max_ps = names(which.max(table(vcf_tab$PS)))
            max_2block_total_coverage = sum(tail(sort(table(vcf_tab$PS)),2)) / n_phased

            if(is.na(block_ratio_thresh)) block_ratio_thresh = max(vcf_tab$block_ratio)
            if(is.na(which_blocks[1])) vcf_tab = vcf_tab[block_ratio >= block_ratio_thresh]
            if(!is.na(which_blocks[1]) & length(which_blocks) > 1) vcf_tab = vcf_tab[block_index %in% which_blocks]

            # max_block_total_coverage = nrow(vcf_tab[block_ratio==max(vcf_tab$block_ratio)]) / n_total_input_SNPs
            max_block_phased_coverage = nrow(vcf_tab[block_ratio==max(vcf_tab$block_ratio)]) / n_phased
            selected_block_phased_coverage = nrow(vcf_tab) / n_phased

            cat('  >> n_total_phased in n_total_input_SNPs:', n_phased/n_total_input_SNPs, ' | n_max_block in n_phased_SNPs:', max_block_phased_coverage, '| selected_blocks in n_phased_SNPs:', selected_block_phased_coverage, '\n')

        }

        if(!is.null(vcf_tab$block_index)) vcf_tab$hap_val = unname(hap_val[vcf_tab$hap])/(vcf_tab$block_index)
        if(is.null(vcf_tab$block_index)) vcf_tab$hap_val = unname(hap_val[vcf_tab$hap])

        if(save2wig==FALSE) return(vcf_tab)

        ## save as wig in the following
        wig_tab = data.table(chr=gsub('chrchr', 'chr', paste0('chr', vcf_tab[[1]])), pos_start=as.character(vcf_tab$POS - 1), pos_end=as.character(vcf_tab$POS), hap=vcf_tab$hap_val)
        wig_tab = na.omit(wig_tab)
        data.table::fwrite(wig_tab, file=wig_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

        if(!phased) return(vcf_tab)

        phased_ground_truth = consist_ratio = NA
        if(length(wig_ground_truth) > 1)
        {
            colnames(wig_ground_truth)[4] = 'hap_ground_truth'
            wig_overlap = merge(wig_tab, wig_ground_truth, by=c('chr', 'pos_start', 'pos_end'))

            phased_ground_truth = nrow(wig_overlap) / nrow(wig_ground_truth)

            snp_coverage = sum(vcf_tab_input$POS %in% wig_ground_truth$pos_end) / nrow(wig_ground_truth)
            phased_coverage = nrow(wig_overlap) / nrow(wig_ground_truth)

            consist_ratio_raw = mean(wig_overlap$hap==wig_overlap$hap_ground_truth)
            consist_ratio = max(consist_ratio_raw, 1 - consist_ratio_raw)
            cat('  >> ground truth snp coverage:', snp_coverage, '\n')
            cat('  >> ground truth phased coverage:', phased_coverage, ' | consist_ratio:', consist_ratio, ' | ')
            consist_sign = sign(consist_ratio_raw - 0.5)
            wig_overlap$hap = consist_sign*(wig_overlap$hap)
            wig_tab_diff = wig_overlap[hap!=hap_ground_truth]
            data.table::fwrite(wig_tab_diff, file=wig_diff_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
        }

        cat('done\n')
        if(consider_ps==TRUE) phasing_info = list(phased_ground_truth=phased_ground_truth, phased_in_total=n_phased/n_total_input_SNPs, max_block_ratio=max_block_phased_coverage, max_2block_ratio=max_2block_total_coverage, consist_ratio=consist_ratio, vcf_tab=vcf_tab)
        if(consider_ps==FALSE) phasing_info = list(phased_ground_truth=phased_ground_truth, max_block_ratio=NA, max_2block_ratio=max_2block_ratio, consist_ratio=consist_ratio, vcf_tab=vcf_tab)
        if(vcf_tab_only) phasing_info = vcf_tab
        return(phasing_info)
    }

    vcf2tab_block = function(vcf_file, gz=FALSE, phased=TRUE, save2wig=FALSE, consider_ps=FALSE, QUAL_thresh=-1, wig_ground_truth=NA, format_shapeit2_hap=FALSE)
    {
        library(data.table)

        if(1)
        {
            if(gz==TRUE)
            {
                wig_file = gsub('.vcf.gz', '.wig', vcf_file, ignore.case=TRUE)
                my_skip = as.numeric(strsplit(system(intern=TRUE, sprintf('zgrep -n "#CHROM\tPOS" %s', vcf_file)), ':')[[1]][1]) - 1
            }

            if(gz==FALSE)
            {
                wig_file = gsub('[.]vcf', '.wig', vcf_file, ignore.case=TRUE)
                my_skip = as.numeric(strsplit(system(intern=TRUE, sprintf('grep -n "#CHROM\tPOS" %s', vcf_file)), ':')[[1]][1]) - 1
            }

            wig_diff_file = gsub('.vcf', '.diff.wig', vcf_file, ignore.case=TRUE)

            vcf_tab = data.table::fread(vcf_file, sep='\t', skip=my_skip)
            vcf_tab = vcf_tab[, mget(setdiff(colnames(vcf_tab), 'INFO'))]
            vcf_tab = vcf_tab[QUAL > QUAL_thresh]
        }

        sample_name = tail(colnames(vcf_tab), 1)
        vcf_tab$hap = substr(vcf_tab[[sample_name]], 1, 3)
        vcf_tab = vcf_tab[hap %in% c('0/1', '0|1', '1|0')[-1]]

        ## is wig_ground_truth is not NA
        if(length(wig_ground_truth) > 1) vcf_tab = vcf_tab[POS %in% wig_ground_truth$pos_end]

        colnames(vcf_tab)[c(1,2)] = c('chr', 'POS') 

        vcf_tab_input = vcf_tab
        n_total_input_SNPs =  nrow(vcf_tab_input) ## total number of variants in the VCF file input for phasing

        #########################################################
        consider_ps = 1
        if(consider_ps) 
        {
            vcf_tab$PS = unname(sapply(vcf_tab[[sample_name]], function(v) tail(strsplit(v, ':')[[1]],1)))
            vcf_tab = vcf_tab[PS!='.']

            n_phased = nrow(vcf_tab)
            max_ps = names(which.max(table(vcf_tab$PS)))
            vcf_tab = vcf_tab[PS==max_ps]

            max_block_total_coverage = nrow(vcf_tab) / n_total_input_SNPs
            max_block_phased_coverage = nrow(vcf_tab) / n_phased

            cat('  >> n_max_block in n_total_input_SNPs:', max_block_total_coverage, ' | n_max_block in n_phased_SNPs:', max_block_phased_coverage, '\n')

        }

        hap_val = c(1,-1, 1)
        names(hap_val) = c('0|1', '1|0', '0/1')
        vcf_tab$hap_val = unname(hap_val[vcf_tab$hap])

        if(save2wig==FALSE) return(vcf_tab)

        if(grepl('chr', vcf_tab[[1]])[1]) wig_tab = data.table(chr=vcf_tab[[1]], pos_start=vcf_tab$POS - 1, pos_end=vcf_tab$POS, hap=vcf_tab$hap)
        if(!grepl('chr', vcf_tab[[1]])[1]) wig_tab = data.table(chr=paste0('chr', vcf_tab[[1]]), pos_start=vcf_tab$POS - 1, pos_end=vcf_tab$POS, hap=vcf_tab$hap)

        wig_tab$hap = hap_val[wig_tab$hap]
        wig_tab = na.omit(wig_tab)

        data.table::fwrite(wig_tab, file=wig_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

        if(!phased) return(vcf_tab)


        if(length(wig_ground_truth) > 1)
        {
            colnames(wig_ground_truth)[4] = 'hap_ground_truth'
            wig_overlap = merge(wig_tab, wig_ground_truth, by=c('chr', 'pos_start', 'pos_end'))

            phased_ground_truth = nrow(wig_overlap) / nrow(wig_ground_truth)

            snp_coverage = sum(vcf_tab_input$POS %in% wig_ground_truth$pos_end) / nrow(wig_ground_truth)
            phased_coverage = nrow(wig_overlap) / nrow(wig_ground_truth)

            consist_ratio_raw = mean(wig_overlap$hap==wig_overlap$hap_ground_truth)
            consist_ratio = max(consist_ratio_raw, 1 - consist_ratio_raw)
            cat('  >> ground truth snp coverage:', snp_coverage, '\n')
            cat('  >> ground truth phased coverage:', phased_coverage, ' | consist_ratio:', consist_ratio, ' | ')
            consist_sign = sign(consist_ratio_raw - 0.5)
            wig_overlap$hap = consist_sign*(wig_overlap$hap)
            wig_tab_diff = wig_overlap[hap!=hap_ground_truth]
            data.table::fwrite(wig_tab_diff, file=wig_diff_file, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
        }

        cat('done\n')
        if(consider_ps==TRUE) phasing_info = list(phased_ground_truth=phased_ground_truth, phased_in_total=n_phased/n_total_input_SNPs, max_block_ratio=max_block_phased_coverage, consist_ratio=consist_ratio, vcf_tab=vcf_tab)
        if(consider_ps==FALSE) phasing_info = list(phased_ground_truth=phased_ground_truth, max_block_ratio=NA, consist_ratio=consist_ratio, vcf_tab=vcf_tab)

        return(phasing_info)
    }


    vcf_filter = function(vcf_raw_file, QUAL_thresh, name_only=FALSE, chr_name=NULL, heto_only=FALSE)
    {
        vcf_filtered_het = paste0(gsub('[.]vcf', '', vcf_raw_file), '.HET.bi.QUAL=', QUAL_thresh, '.vcf')
        vcf_filtered_all = paste0(gsub('[.]vcf', '', vcf_raw_file), '.ALL.bi.QUAL=', QUAL_thresh, '.vcf')

        if(name_only) return(vcf_filtered_het) ## only return the name
 
        command_filter_het = sprintf('%s view --max-alleles 2 %s | %s filter -i \'GT="het" && QUAL > %s\' > %s ', bcftools, vcf_raw_file, bcftools, QUAL_thresh, vcf_filtered_het)
        command_filter_all = sprintf('%s view --max-alleles 2 %s | %s filter -i \'QUAL > %s\' > %s ', bcftools, vcf_raw_file, bcftools, QUAL_thresh, vcf_filtered_all)

        system(command_filter_het)
        if(!heto_only) system(command_filter_all)

        count_raw = as.numeric(strsplit(system(intern=TRUE, sprintf('wc -l %s', vcf_raw_file)), ' ')[[1]][[1]])
        count_filtered = as.numeric(strsplit(system(intern=TRUE, sprintf('wc -l %s', vcf_filtered_het)), ' ')[[1]][[1]])
        cat('count_raw:', count_raw, ' | count_filtered_het:', count_filtered, '\n')
        return(vcf_filtered_het)
    }
    
    ######################################################### VCF by freebayes

    stat_ins = function(wk_dir, no_ins=FALSE, chrs) ## stat of over-represented sequences ## insertion
    {
        chr_chars = paste0('chr', chrs)
        n_snp_plot_file = file.path(wk_dir, 'mega/VCFs/', 'variants_type_HET.pdf')
        if(no_ins) n_snp_plot_file = file.path(wk_dir, 'mega/VCFs/', 'variants_type_HET.no_ins.pdf')
        stat_tab_ALL = foreach(chr_char=chr_chars, .combine=rbind) %dopar%
        {
            vcf_raw_file = sprintf("%s/%s/mega.%s.HET.bi.QUAL=10.vcf", wk_dir, 'mega/VCFs/', chr_char)
            command_stat = sprintf("%s view %s -H | cut -f5 | cut -c 2- | sort | uniq -c | sort", bcftools, vcf_raw_file)
            stat = trimws( system(command_stat, intern=TRUE) )

            if(length(stat) ==1) {stat_tab = data.table(n_count=as.numeric(stat), snp='SNV', chr=chr_char); return(stat_tab)}
            if(length(stat) ==0) {stat_tab = data.table(n_count=0, snp='SNV', chr=chr_char); return(stat_tab)}

            stat_tab = data.table( do.call(rbind, strsplit(stat, " ")) )
            colnames(stat_tab) = c("n_count", 'snp')
            stat_tab[n_count==snp]$snp = "SNV" 
            stat_tab$chr = chr_char
            stat_tab$n_count = as.numeric(stat_tab$n_count)
            stat_tab
            # stat_tab[nchar(snp) > 2]
        }

        snp_count = tapply(stat_tab_ALL$n_count, stat_tab_ALL$snp, sum)
        snp_count_tab = data.table(snp=names(snp_count), n_count=unname(snp_count))
        snp_count_tab$ratio = snp_count_tab$n_count / sum(snp_count_tab$n_count )
        snp_count_tab = snp_count_tab[order(-n_count)][1:10]
        snp_count_tab$rank = 1:nrow(snp_count_tab)

        snp2rm = snp_count_tab[ratio > 0.01 & nchar(snp) > 2 & snp!="SNV" ]$snp

        pdf(n_snp_plot_file, width=8, height=6)
        p = ggplot(data=snp_count_tab, aes(x=rank, y=ratio, label=substr(ratio, 1, 6))) + geom_line(color="red") + geom_point()
        p = p + geom_label() + scale_x_continuous("", labels=snp_count_tab$snp, breaks=1:10)
        print(p)    
        graphics.off()

        return(snp2rm)
    }

    parallel_free_v2 = function(free_exe, ref_fa, bam_head, wk_dir, chrs)
    {
        vcf_prefix = gsub('pos_sorted', '', bam_head)
        QUAL_thresh = 10
        
        no_out = foreach(chr=chrs) %dopar%
        {
            pos_sorted_bam_chr = sprintf("%s.chr%s.bam", bam_head, chr)
            vcf_des_file = gsub('splits', 'VCFs', paste0(vcf_prefix, 'chr', chr, '.vcf'))
            # command_free =  sprintf('%s -f %s %s -r chr%s | vcfallelicprimitives -kg > %s', free_exe, ref_fa, pos_sorted_bam_chr, chr, vcf_des_file)
            command_free =  sprintf('%s -f %s %s -r chr%s -m 10 --min-coverage 5 -g 1000 | vcfallelicprimitives -kg > %s', free_exe, ref_fa, pos_sorted_bam_chr, chr, vcf_des_file)
            cat(command_free, '\n')
            system(command_free)
            vcf_filter(vcf_raw_file=vcf_des_file, QUAL_thresh=QUAL_thresh, name_only=FALSE)
            # unlink(vcf_des_file)
        }

        snp2rm = stat_ins(wk_dir, no_ins=FALSE, chrs) ## insertion
        snp2rm = c(snp2rm, 'GATC') ## add GATC, as it freaquently appears in MboI and Arima

        if(length(snp2rm) <= 0) return()
        snp2rm_pattern = paste0(snp2rm, collapse="\\|")

        ############################### else remove those lines containing such snp2rm

        no_out = foreach(chr=chrs) %dopar% ## do it again but only to keep snv and remove ins
        {
            vcf_des_file = gsub('splits', 'VCFs', paste0(vcf_prefix, 'chr', chr, '.vcf'))
            vcf_filtered_het = paste0(gsub('[.]vcf', '', vcf_des_file), '.HET.bi.QUAL=', QUAL_thresh, '.vcf')
            vcf_filtered_all = paste0(gsub('[.]vcf', '', vcf_des_file), '.ALL.bi.QUAL=', QUAL_thresh, '.vcf') 
            command_filter_het = sprintf("sed -i '/%s/d' %s", snp2rm_pattern, vcf_filtered_het)
            command_filter_all = sprintf("sed -i '/%s/d' %s", snp2rm_pattern, vcf_filtered_all)
            system(command_filter_het)
            system(command_filter_all)
        }

        try(silent=TRUE, stat_ins(wk_dir, no_ins=TRUE, chrs))
    }

    run_juicer = function(wk_dir, run_name, fastq_prefix, STEP, thread, juicerDir, enzyme, res_sites, ref_fa, chrom_sizes, server_capacity_low=FALSE, no_frag=0)
    {
        # system('export PATH=/mnt/ed4/Yuanlong/3.Packages/bwa/:$PATH')
        res_dir = file.path(wk_dir, run_name)
        logfile = file.path(wk_dir, run_name, 'log.txt')
        time0 = Sys.time()

        #################################

        if(grepl('mm10', ref_fa)) genome_version = 'mm10'
        if(grepl('hg19', ref_fa)) genome_version = 'hg19'

        cmd = sprintf("%s/scripts/my_juicer.sh -k %s -s %s -n %s -g %s -z %s -d %s -D %s -y %s -p %s -t %s >> %s", juicerDir, STEP, enzyme, fastq_prefix, genome_version, ref_fa, res_dir, juicerDir, res_sites, chrom_sizes, thread, logfile)
        if(server_capacity_low) cmd = sprintf("%s/scripts/my_juicer_gravitron.sh -k %s -s %s -n %s -g %s -z %s -d %s -D %s -y %s -p %s -t %s >> %s", juicerDir, STEP, enzyme, fastq_prefix, genome_version, ref_fa, res_dir, juicerDir, res_sites, chrom_sizes, thread, logfile)

        if(no_frag==1) cmd = sprintf("%s/scripts/my_juicer_no_frag.sh -k %s -s %s -n %s -g %s -z %s -d %s -D %s -y %s -p %s -t %s >> %s", juicerDir, STEP, enzyme, fastq_prefix, genome_version, ref_fa, res_dir, juicerDir, res_sites, chrom_sizes, thread, logfile)

        cat(cmd)
        cat(as.character(Sys.time()), '\n', file=logfile)
        cat('cmd:', cmd, file=logfile, append=TRUE)  
        system(cmd)

        if(STEP==5)
        {
            hic_file_ori = sprintf("%s/%s/aligned/inter.hic", wk_dir, run_name)
            hic_file_rename = sprintf("%s/%s/aligned/%s.hic", wk_dir, run_name, gsub("-", "_", fastq_prefix))
            try(silent=TRUE, system(sprintf("mv %s %s", hic_file_ori, hic_file_rename)))
        }        
    }

    run_juicer_mega = function(wk_dir, juicerDir, enzyme, res_sites)
    {
        if(grepl('mm10', res_sites)) genome_version = 'mm10'
        if(grepl('hg19', res_sites)) genome_version = 'hg19'
        cmd = sprintf("%s/scripts/my_mega.sh -s %s -y %s -g %s -D %s -d %s", juicerDir, enzyme, res_sites, genome_version, juicerDir, wk_dir)
        system(cmd)
    }


    vcf_filter_by_QUAL = function(chr_char, QUAL_thresh)
    {
        vcf_chr_raw = sprintf('%s/%s.%s.vcf', VCF_DIR, VCF_prefix, chr_char)
        vcf_chr = vcf_filter(vcf_chr_raw, QUAL_thresh=QUAL_thresh, name_only=FALSE) ## this line should not be included to the full dopar loop, with be overwritten
    }

    ##############################

    pre = function(mat2pre, hic_file, ref_ref_fa="hg19")
    {
        file2pre = paste0(hic_file, ".delete")
        fwrite(mat2pre, file=file2pre, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
        command_pre = paste0('java -Xmx200g -jar ', juicer_jar, ' pre -r 10000,20000,25000,40000,50000,100000,500000 ', file2pre, ' ', hic_file, " ", ref_ref_fa, ' -d true')
        system(command_pre)
        unlink(file2pre)
    }


    ####################################### generate chr_pos and mat_pat for shapeit

    generate_chr_pos_mat_pat_shapeit = function(wk_dir, chr_num)
    {
        chr_name = paste0("chr", chr_num)
        haps = fread(sprintf("%s/mega/HapCut/shapeit_out/mega.chr%s.HET.bi.QUAL=10.vcf.haps", wk_dir, chr_num))
        haps = haps[order(V3)]
        chr_pos_shapeit = sprintf("%s/mega/HapCut/phased_hap_intg/hic_chr%s_shapeit.chr_pos", wk_dir, chr_num)
        pat_mat_shapeit = sprintf("%s/mega/HapCut/phased_hap_intg/hic_chr%s_shapeit.pat_mat", wk_dir, chr_num)

        chr_pos = data.table(cbind(chr_name, t(haps[,3])))
        pat_mat = data.table(chr=paste0(chr_name, ":", haps[[3]]), mat=ifelse(haps$V7==1, haps$V4, haps$V5), pat=ifelse(haps$V7==1, haps$V5, haps$V4))

        fwrite(chr_pos, file=chr_pos_shapeit, row.names=FALSE, col.names=FALSE, sep=" ")
        fwrite(pat_mat, file=pat_mat_shapeit, row.names=FALSE, col.names=FALSE, sep=" ")
    }


    merge_bams = function(bams_in, mega_bam) ## sorted by name
    {
        dir.create(dirname(mega_bam), showWarnings=FALSE, recursive=TRUE)
        tmp_dir = '/scratch/yliu5/4.SRA_tmp/'
        command_merge = sprintf('samtools merge -f -@ 30 -n %s %s', paste0(mega_bam, '.tmp'), paste0(bams_in, collapse=" "))
        cat(command_merge, file=gsub(".bam", ".merge.log", mega_bam))
        
        ## this memory setting should be different for graviton and positron
        command_merge_sort = sprintf('samtools sort -@ 15 -m 2G %s -O BAM -o %s -T %s/tmp', paste0(mega_bam, '.tmp'), mega_bam, tmp_dir)
        command_index =  sprintf('samtools index -@ 20 %s', mega_bam)

        system(command_merge)
        system(command_merge_sort)
        unlink(paste0(mega_bam, '.tmp'))
        system(command_index)
    }

    merge_bams_sorted = function(bams_in, mega_bam) ## sorted by coordinates
    {
        dir.create(dirname(mega_bam), showWarnings=FALSE, recursive=TRUE)
        tmp_dir = '/scratch/yliu5/4.SRA_tmp/'
        command_merge = sprintf('samtools merge -f -@ 30 %s %s', paste0(mega_bam, '.tmp'), paste0(bams_in, collapse=" "))
        cat(command_merge, file=gsub(".bam", ".merge.log", mega_bam))
        
        ## this memory setting should be different for graviton and positron
        command_merge_sort = sprintf('samtools sort -@ 15 -m 2G %s -O BAM -o %s -T %s/tmp', paste0(mega_bam, '.tmp'), mega_bam, tmp_dir)
        command_index =  sprintf('samtools index -@ 20 %s', mega_bam)

        system(command_merge)
        system(command_merge_sort)
        unlink(paste0(mega_bam, '.tmp'))
        system(command_index)
    }

    plot_colorByDensity = function(x1,x2,
                               ylim=c(min(x2),max(x2)),
                               xlim=c(min(x1),max(x1)),
                               xlab="", ylab="", main="", cex=2, ...) {
     
    df <- data.frame(x1,x2)
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
    # cols = bluered(1000)
    df$col <- cols[df$dens]
    plot(x2~x1, data=df[order(df$dens),], 
         ylim=ylim,xlim=xlim,pch=20,col=col,
         cex=cex,xlab=xlab,ylab=ylab,
         main=main, ...)
    }

   
    split_fastq_helper = function(fastq_prefix, n_part)
    { 
        wk_dir = dirname(fastq_prefix)
        base = basename(fastq_prefix)

        dir.create(file.path(wk_dir, 'fastq_bp'), recursive=TRUE, showWarnings=FALSE)
        cmd_split = sprintf( "%s split2 -j 5 -1 %s_1.fastq.gz -2 %s_2.fastq.gz -p %s -O %s -f", seqkit, fastq_prefix, fastq_prefix, n_part, fastq_prefix )
        out_flag = system(cmd_split)
        cat('split done\n')

        no_out = foreach(part=1:n_part) %do%
        {
            no_out = foreach(mate=1:2) %do%
            {
                id = sprintf('00%s', part)
                id = substr(id, nchar(id) - 2, nchar(id))
                f_old = sprintf("%s/%s_%s.part_%s.fastq.gz", fastq_prefix, base, mate, id)
                f_new = sprintf("%s/%s.PART%s_%s.fastq.gz", wk_dir, base, part, mate) 
                file.rename(from=f_old, to=f_new)             
            }
        }

        unlink(fastq_prefix, recursive=TRUE) ## just remove the split folder

        no_out=foreach(mate=1:2) %do% ## backup the original file before deciding to remove
        {
            f_old = sprintf("%s_%s.fastq.gz", fastq_prefix, mate)
            f_new = sprintf("%s/fastq_bp/%s_%s.fastq.gz", wk_dir, base, mate)
            file.rename(from=f_old, to=f_new)
        }

        return(out_flag)
    }


    split_fastq = function(wk_dir, fastq_id, sizeGb=5)
    {
        fastq_prefix = file.path(wk_dir, fastq_id)
        f1 = paste0(fastq_prefix, '_1.fastq.gz')
        f2 = paste0(fastq_prefix, '_2.fastq.gz')

        if(!file.exists(f2)) return(NULL)
        f_sizes = as.numeric(do.call(rbind, strsplit(utils:::format.object_size(file.size(f1), "Gb"), ' '))[,1])
        n_parts = max( ceiling( f_sizes / sizeGb ) ) ## split if file greater than 8Gb
        if(n_parts < 2) 
        {
            cat("fastq file is already small, no need to split\n")
            return(NULL)
        }

        split_fastq_helper(fastq_prefix, n_part=n_parts)
    }


    split_fastq_helper_dir = function(wk_dir) ## split all files in the folder
    {
        fs = list.files(wk_dir, 'fastq.gz', full.names=TRUE)
        f_sizes = as.numeric(do.call(rbind, strsplit(utils:::format.object_size(file.size(fs), "Gb"), ' '))[,1])
        n_parts = ceiling( f_sizes / 8 ) ## split if file greater than 8Gb
        split_tab = data.table(f=fs, n_parts=n_parts, f_size=f_sizes)[n_parts > 1]
        split_tab$fastq_prefix = gsub('_2.fastq.gz', '', gsub('_1.fastq.gz', '', split_tab$f))
        split_tab = unique(split_tab, by=c('fastq_prefix'))

        fastq_prefix2split = unique( split_tab$fastq_prefix )

        ##############################

        registerDoParallel(cores = length(fastq_prefix2split))
        no_out = foreach(i=1:nrow(split_tab)) %dopar% 
        {
            cat(sprintf('Begin split: %s', split_tab$fastq_prefix[i]), '\n')
            split_fastq_helper(split_tab$fastq_prefix[i], n_part=split_tab$n_parts[i])
            cat(sprintf('End split: %s', split_tab$fastq_prefix[i]), '\n')
        }
    }


    # get_chr2use = function(wk_dir) ## do not analyze chrs with less than 1E3 SNPs
    # {
    #     chr2phase_f = file.path(wk_dir, 'mega/VCFs/chr2phase.txt')
    #     if(!file.exists(chr2phase_f))
    #     {
    #         snp_count_f = file.path(wk_dir, 'mega/VCFs/n_variants_het.no_ins.txt')
    #         snp_count_tab = fread(snp_count_f)
    #         colnames(snp_count_tab) = c('chr', 'n', 'QUAL')
    #         chr_num2use = as.character(unique(snp_count_tab[QUAL==1][n > 1E3]$chr)) ## IMR90 chr22 is 5E4
    #         chr_num2use[chr_num2use=='23'] = 'X'
    #         cat(chr_num2use, file=chr2phase_f)
    #     }
    #     chr_num2use = unname(unlist(fread(chr2phase_f)))
    #     return(chr_num2use)
    # }

    #########################################################
    
    # my_forceSymmetric = function(mat) ## Matrix package does not work on UNIL HPC conda env intg
    # {
    #     new_mat = mat + t(mat)
    #     diag(new_mat) = diag(mat)
    #     return(new_mat)
    # }
 
    #########################################################

    # get_fragile_sites = function(vcf_tab_li, save_file) ## I will not segment further inside blocks. This can be further optimized
    # {
    #     my_merge = function(...) merge(..., by='POS', all=TRUE)
    #     vcf_tab_li_slim = vcf_tab_li
    #     len = length(vcf_tab_li_slim)
    #     n_col = len + 1

    #     for(i in 1:len)
    #     {
    #         n_blocks = length(unique(vcf_tab_li[[i]]$block_index))
    #         print(n_blocks)

    #         set.seed(i) ## for different vcfs
    #         rand_vals = round(abs(rnorm(n_blocks, 10, 10)), 5) ## add a random value for the same block
    #         vcf_tab_li_slim[[i]]$hap_val = (vcf_tab_li_slim[[i]]$hap_val)*(rand_vals[(vcf_tab_li_slim[[i]]$block_index)]) ## sam abs value if coming from the same hap block
    #         vcf_tab_li_slim[[i]] = vcf_tab_li_slim[[i]][, mget(c('POS', 'hap_val'))]
    #         colnames(vcf_tab_li_slim[[i]])[2] = paste0(colnames(vcf_tab_li_slim[[i]])[2], '_', names(vcf_tab_li)[i])
    #     }

    #     vcf_combined = Reduce(my_merge, vcf_tab_li_slim)
    #     vcf_combined = vcf_combined[order(POS)]

    #     # signs = sign(vcf_combined$hap_val_p1_1000)

    #     na_count = apply(vcf_combined[,-1], 2, function(v) sum(is.na(v)) )
    #     low_na_para = names(which.min(na_count))

    #     signs = unname(unlist(sign(vcf_combined[, mget(low_na_para)])))
    #     for(i in 2:ncol(vcf_combined)) vcf_combined[[i]] = signs*vcf_combined[[i]] ## set sign of the longest phased as +1  

    #     vcf_combined$ALL =  apply(abs(vcf_combined[, -c(1, 2), with=FALSE] - vcf_combined[, -c(1, n_col), with = FALSE]), 1, mean, na.rm=TRUE) ## combinations
    #     vcf_combined = vcf_combined[!is.nan(ALL)]
    #     vcf_combined[ALL<1]$ALL = vcf_combined[ALL<1]$ALL + 1

    #     vcf_combined[, `:=` (segment_ratio = .N/nrow(vcf_combined)), by = ALL] ## ratio of each segment
    #     # vcf_combined = vcf_combined[segment_ratio > 0.05]
    #     vcf_combined = vcf_combined[segment_ratio > 0.02] ## changed to 0.01 on 03 Sep 2023; better to have smaller segments
    #     vcf_combined$segment_ratio = substr(vcf_combined$segment_ratio, 1, 10)

    #     centoids_raw = tapply(vcf_combined$POS/1E6, vcf_combined$segment_ratio, median) ## center pos
    #     n_segments = length(centoids_raw)
    #     segments = 1:n_segments
    #     names(segments) = names(sort(centoids_raw))
    #     vcf_combined$segment_index = unname(segments[vcf_combined$segment_ratio])
    #     centoids = tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, median)
    #     # fivenums = do.call(rbind, tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, summary))
    
    #     ten_percentile = c(sort(tapply(vcf_combined$POS/1E6, vcf_combined$segment_index, quantile, 0.01)*1E6), Inf)
    #     uniq_ALLs = unique(vcf_combined$ALL) ## already ordered based on first appearance of POS


    #     vcf_combined$stair_val = 0
    #     for(z in 1:length(uniq_ALLs)) vcf_combined[between(POS, ten_percentile[z], ten_percentile[z+1])]$stair_val = apply(vcf_combined[between(POS, ten_percentile[z], ten_percentile[z+1]), 'ALL'], 1, function(v) sum(setdiff(uniq_ALLs[1:z], v)))
        
    #     # vcf_combined$stair_val = substr(vcf_combined$stair_val, 1, 20)
    #     vcf_combined[, `:=` (new_count = .N/nrow(vcf_combined)), by = stair_val]
    #     # vcf_combined = vcf_combined[new_count > 0.05] ## ignore those with less than 5%
    #     vcf_combined = vcf_combined[new_count > 0.02] ## ignore those with less than 5%        
    #     vcf_combined$segment_index = match(vcf_combined$stair_val, unique(vcf_combined$stair_val))

    #     n_segments = length(unique(vcf_combined$segment_index))
    #     fragile_sites = NULL
    #     for(s in 1:(n_segments-1)) 
    #     {
    #         quant_1 = quantile(vcf_combined[segment_index==s]$POS, 0.95) 
    #         quant_2 = quantile(vcf_combined[segment_index==(s+1)]$POS, 0.05)
    #         mid = (quant_1 + quant_2)/2/1E6
    #         fragile_sites = c(fragile_sites, mid)
    #     }

    #     fragile_info = list(fragile_sites=fragile_sites, vcf_combined=vcf_combined)

    #     pdf(save_file, width=10, height=6)
    #     plot(xaxt="n", xlab='Mb', ylab='segment pseudo val', vcf_combined$POS/1E6, vcf_combined$ALL, pch='.', col='red', cex=0.1, type='p')
    #     axis(1, at=(1:100)*5, las=2)
    #     abline(v=fragile_sites, lty=2, col='blue')
    #     graphics.off()

    #     return(fragile_info)
    # }


  ## candidate 1 / 13Nov2023
   # get_fragile_sites = function(vcf_tab_li, save_file) ## I will not seg further inside blocks. This can be further optimized
   #  {
   #      my_merge = function(...) merge(..., by='POS', all=TRUE)
   #      vcf_tab_li_slim = vcf_tab_li
   #      len = length(vcf_tab_li_slim)
   #      n_col = len + 1

   #      for(i in 1:len)
   #      {
   #          n_blocks = length(unique(vcf_tab_li[[i]]$block_index))
   #          print(n_blocks)

   #          set.seed(i) ## for different vcfs
   #          rand_vals = round(abs(rnorm(n_blocks, 10, 10)), 5) ## add a random value for the same block
   #          vcf_tab_li_slim[[i]]$hap_val = (vcf_tab_li_slim[[i]]$hap_val)*(rand_vals[(vcf_tab_li_slim[[i]]$block_index)]) ## sam abs value if coming from the same hap block
   #          vcf_tab_li_slim[[i]] = vcf_tab_li_slim[[i]][, mget(c('POS', 'hap_val'))]
   #          colnames(vcf_tab_li_slim[[i]])[2] = paste0(colnames(vcf_tab_li_slim[[i]])[2], '_', names(vcf_tab_li)[i])
   #      }

   #      vcf_merge = Reduce(my_merge, vcf_tab_li_slim)
   #      vcf_merge = vcf_merge[order(POS)]

   #      na_count = apply(vcf_merge[,-1], 2, function(v) sum(is.na(v)) )
   #      low_na_para = names(which.min(na_count))

   #      signs = unname(unlist(sign(vcf_merge[, mget(low_na_para)])))
   #      for(i in 2:ncol(vcf_merge)) vcf_merge[[i]] = signs*vcf_merge[[i]] ## set sign of the longest phased as +1  

   #      vcf_merge$ALL =  apply(abs(vcf_merge[, -c(1, 2), with=FALSE] - vcf_merge[, -c(1, n_col), with = FALSE]), 1, mean, na.rm=TRUE) ## combinations
   #      vcf_merge = vcf_merge[!is.nan(ALL)]
   #      vcf_merge[ALL<1]$ALL = vcf_merge[ALL<1]$ALL + 1

   #      vcf_merge[, `:=` (group_ratio = .N/nrow(vcf_merge)), by = ALL] ## ratio of each seg
   #      vcf_merge = vcf_merge[group_ratio > 0.05] ## changed to 0.01 on 03 Sep 2023; better to have smaller groups
   #      vcf_merge$group_ratio = substr(vcf_merge$group_ratio, 1, 10)
    
   #      #####################################

   #      n_groups = length(unique(vcf_merge$ALL))
   #      vcf_merge$group_index = match(as.numeric(vcf_merge$group_ratio), sort(decreasing=TRUE, unique(as.numeric(vcf_merge$group_ratio))))
   #      vcf_merge$group_index = factor(vcf_merge$group_index, levels=1:n_groups)

   #      vcf_merge[, group_size:=length(POS), by='group_index']
   #      group_size_tab = unique(vcf_merge, by='group_index')[, c('group_index', 'group_size')]
   #      colnames(group_size_tab)[1] = 'group'

   #      #####################################

   #      p = ggplot(vcf_merge, aes(x=POS/1E6, color=group_index)) + geom_density()
   #      plot_data = as.data.table(ggplot_build(p)$data[[1]])

   #      dt = plot_data[, c('x', 'y', 'group')][order(group, x)]
   #      dt$group = factor(dt$group, levels=1:n_groups)
   #      dt = merge(dt, group_size_tab, by='group')
   #      dt$y = dt$y * dt$group_size ## multiplied by group size

   #      dt_wide = dcast(dt, x ~ group, value.var = "y")
   #      dt_wide$max_group = apply(dt_wide[,-1], 1, which.max)

   #      ##################################### find f sites

   #      density_breaks = which(diff(dt_wide$max_group)!=0)
   #      fragile_sites = dt[group==1]$x[density_breaks]
   #      fragile_sites = c(0, fragile_sites, Inf)
        
   #      # refine_fsites = function(dens_thresh=100) ## remove sites that are too nearby and with low SNP density between, thus to reduce fragile_sites
   #      {  
   #          dens_thresh = 100
   #          k = 1
   #          while(any(abs(diff(fragile_sites)) < 10)) ## if any fragile_sites less than 10Mb
   #          {
   #              n_fsites = length(fragile_sites)
   #              max_density_between_fs = foreach(i=1:(length(fragile_sites)-1), .combine=c) %do%
   #              {
   #                  max(dt[between(x, fragile_sites[i], fragile_sites[i+1])]$y)
   #              }

   #              max_density_at_fs = foreach(i=1:length(fragile_sites), .combine=c) %do%
   #              {
   #                  tmp = dt[between(x, fragile_sites[i], fragile_sites[i])]
   #                  if(nrow(tmp)==0) return(0)
   #                  max(tmp$y)
   #              }

   #              f_tab = data.frame(f_site=fragile_sites[1:(n_fsites-1)], max_between=max_density_between_fs, max_at=max_density_at_fs[1:(n_fsites-1)])
                
   #              if(min(f_tab$max_between) > dens_thresh) break
   #              if(k > 500) break
   #              if(length(fragile_sites) <= 3) break

   #              min_index = which.min(f_tab$max_between)
   #              index2rm_tmp = c(min_index, min_index+1)
   #              index2rm = index2rm_tmp[which.max(f_tab$max_at[index2rm_tmp])]
   #              f_tab = f_tab[-index2rm, ]
   #              fragile_sites = unique(c(0, f_tab$f_site, Inf))
   #              k = k + 1
   #          }
   #      }

   #      # refine_fsites = function(dens_thresh=100) ## remove sites that are too nearby and with low SNP density between, thus to reduce fragile_sites

   #      {  
   #          while(any(abs(diff(fragile_sites)) < 1)) ## if any fragile_sites less than 1Mb
   #          {
   #              n_fsites = length(fragile_sites)

   #              max_density_at_fs = foreach(i=1:length(fragile_sites), .combine=c) %do%
   #              {
   #                  tmp = dt[between(x, fragile_sites[i], fragile_sites[i])]
   #                  if(nrow(tmp)==0) return(0)
   #                  max(tmp$y)
   #              }

   #              f_tab = data.frame(f_site=fragile_sites[1:(n_fsites-1)], dist=c(1E7, abs(diff(fragile_sites)))[1:(n_fsites-1)], max_at=max_density_at_fs[1:(n_fsites-1)])
                
   #              if(k > 500) break
   #              if(length(fragile_sites) <= 3) break

   #              min_index = which.min(f_tab$dist)
   #              index2rm_tmp = c(min_index-1, min_index)
   #              index2rm = index2rm_tmp[which.max(f_tab$max_at[index2rm_tmp])]
   #              f_tab = f_tab[-index2rm, ]
   #              fragile_sites = unique(c(0, f_tab$f_site, Inf))
   #              k = k + 1
   #          }
   #      }


    #     fragile_sites = setdiff(fragile_sites, c(0, Inf))

    #     #####################################

    #     z = ggplot(dt, aes(x=x, y=y)) + geom_line(aes(color=group, alpha=0.5))
    #     z = z + scale_x_continuous(breaks=5*(1:100)) + scale_color_manual(values=viridis(n_groups)) + theme_classic() + theme(legend.position='none', plot.margin = margin(3, 1, 1, 1, "cm"), axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
    #     z = z + xlab('Mb') + ylab('snp count') + geom_vline( xintercept=fragile_sites, linetype='dashed', col='blue')
    #     z = z + coord_cartesian(ylim=c(0, min(max(dt$y), 5000)))

    #     n_groups = length(unique(vcf_merge$group_index))
    #     fragile_info = list(fragile_sites=fragile_sites, vcf_merge=vcf_merge)

    #     pdf(save_file, width=10, height=4)
    #     print(z)
    #     plot_colorByDensity(xaxt="n", xlab='Mb', ylab='seg pseudo val', vcf_merge$POS/1E6, vcf_merge$ALL, cex=0.1)
    #     axis(1, at=(1:100)*5, las=2)
    #     abline(v=fragile_sites, lty=2, col='blue')
    #     graphics.off()

    #     return(fragile_info)
    # }
