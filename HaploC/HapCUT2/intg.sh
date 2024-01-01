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

    mkdir -p $wk_dir/mega/HapCut/fragment_hapcut
    mkdir -p $wk_dir/mega/HapCut/phased_hap_hapcut
    mkdir -p $wk_dir/mega/HapCut/hic_phased

    if [ $genome_version == hg19 ]; then
       mkdir -p $wk_dir/mega/HapCut/phased_hap_intg
       mkdir -p $wk_dir/mega/HapCut/fragment_intg
       mkdir -p $wk_dir/mega/HapCut/shapeit_out
    fi

    ################################################################################

    intg_fun_mini = function(chunk_names, chr_char, QUAL_thresh, maxIS_Mb, skip_if_exist=TRUE)
    {
        frag_hapcut_chr = sprintf('%s/mega/HapCut/fragment_hapcut/hic_%s_QUAL=%s_maxIS=%s.fragment', wk_dir, chr_char, QUAL_thresh, maxIS_Mb)
        frag_intg_chr = sprintf('%s/mega/HapCut/fragment_intg/hic_%s_QUAL=%s_maxIS=%s.fragment', wk_dir, chr_char, QUAL_thresh, maxIS_Mb)

        vcf_chr_raw = sprintf('%s/%s.%s.vcf', VCF_DIR, VCF_prefix, chr_char)
        vcf_chr = vcf_filter(vcf_chr_raw, QUAL_thresh=QUAL_thresh, name_only=TRUE) ## this line should not be included to the full dopar loop, will be overwritten
        vcf_chr = gsub("HET", "ALL", vcf_chr) ## because using ALL instead of HET. Otherwise the program will break I think

        # Construct file paths
        frag_hapcut_chr="${wk_dir}/mega/HapCut/fragment_hapcut/hic_${chr_char}_QUAL=${QUAL_thresh}_maxIS=${maxIS_Mb}.fragment"
        frag_intg_chr="${wk_dir}/mega/HapCut/fragment_intg/hic_${chr_char}_QUAL=${QUAL_thresh}_maxIS=${maxIS_Mb}.fragment"
        vcf_chr_raw="${VCF_DIR}/${VCF_prefix}.${chr_char}.vcf"


        extracthairs_chunks()
        {
            for (( i=0; i<${#chunk_names[@]}; i++ )); do
                chunk_num=$((i+1))
                chunk_name=${chunk_names[$i]}
                bam_chr="${wk_dir}/${chunk_name}/hic_separated/hic_${chr_char}.bam"

                ## Construct and execute the command
                ## --mbq is 13 by default -10*log10(0.05) = 13
                ## HiC or hic does not matter, 23 June 2023

                if [ "$chunk_num" -eq 1 ]; then
                    command_fragment="${EXTRACTHAIRS} --HiC 1 --bam ${bam_chr} --VCF ${vcf_chr} --maxIS $(($maxIS_Mb*1000000)) --indels 1 --mmq 10 --mbq 10 --ref ${ref_fa} > ${frag_hapcut_chr}"
                else
                    command_fragment="${EXTRACTHAIRS} --HiC 1 --bam ${bam_chr} --VCF ${vcf_chr} --maxIS $(($maxIS_Mb*1000000)) --indels 1 --mmq 10 --mbq 10 --ref ${ref_fa} >> ${frag_hapcut_chr}"
                fi
                eval $command_fragment
            done
        }

        # Construct the output file path
        output_hap_intg=${wk_dir}/mega/HapCut/phased_hap_intg/hic_${chr_char}_QUAL=${QUAL_thresh}_maxIS=${maxIS_Mb}.hap

        # Check if file exists and its size
        if [ -f "$output_hap_intg" ]; then
            fz=$(du -m "$output_hap_intg" | cut -f1)  # File size in MB
            if [[ ! -z "$fz" && "$fz" > 0.1 && "$skip_if_exist" = true ]]; then
                echo "File exist: $output_hap_intg"
                exit 0  # Equivalent to 'return(NA)' in R
            fi
        fi

        {

        # Check genome version and execute commands
        if [ $genome_version == hg19 ]; then
            # Construct and execute the SHAPEIT2 command
            command_shapeit2="$SHAPEIT2 $chr_num $vcf_chr $frag_hapcut_chr $frag_intg_chr $shapeit_out_dir append"
            eval $command_shapeit2

            # Construct and execute the HapCUT22 command
            output_hap_intg=$wk_dir/mega/HapCut/phased_hap_intg/hic_${chr_char}_QUAL=${QUAL_thresh}_maxIS=${maxIS_Mb}.hap
            output_model=$wk_dir/mega/HapCut/phased_hap_intg/hic_${chr_char}_QUAL=${QUAL_thresh}_maxIS=${maxIS_Mb}.htrans_model
            command_intg="$HapCUT22 --fragments $frag_intg_chr --vcf $vcf_chr --output $output_hap_intg --hic 1 --htrans_data_outfile $output_model --outvcf 1"
            eval $command_intg
        fi
        }
    }

    ################################################################################ intg

    no_out = foreach(k=1:nrow(para_tab)) %dopar%
    {
        chr_char = para_tab$chr_char[k]
        QUAL_thresh = para_tab$QUAL_thresh[k]
        maxIS_Mb = para_tab$maxIS_Mb[k]

        try(intg_fun_mini(chunk_names=chunk_names, chr_char, QUAL_thresh, maxIS_Mb, skip_if_exist=FALSE))
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
