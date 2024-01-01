#!/usr/bin/awk -f    
##########
#The MIT License (MIT)
# Header of VCF
BEGIN {
    print "No genotype information available";

    OFS="\t";
}
$1 ~ /^#CHROM/ {
    # header line
    if (NF <= 8) {
        print "No genotype information available";
        exit;
    }
    else {
        for (i=10; i<=NF; i++){
            samples[i]=$i;
        }
    }
    prev="notset";
}
$0 !~ /^#/ {
    # grab the sample fields and parse phasing
    for (i=10; i<=NF; i++){
        split($i, a, ":");
        split($5, alt, ",");
        toolong=0;
        for (j=1; j <=length(alt); j++) {
            if (length(alt[j]) > 1) toolong=1;
        }
        # only interested in phased SNPs, not indels
        if (a[1] ~ /\|/ && length($4)==1 && !toolong) {
            # position 0 in phase corresponds to ref
            alt[0]=$4;
            split(a[1], gt, "|");
            # phase is different
            if (gt[1] != gt[2]) {
                # first allele is "paternal" and listed first
                print $1":"$2, alt[gt[1]], alt[gt[2]] >> samples[i]"_paternal_maternal.txt";
                # new chromosome
                if ($1 != prev) {
                    if (prev =="notset") {
                        printf("%s ", $1) >> samples[i]"_chr_pos.txt";
                    }
                    else {
                        printf("\n%s ", $1) >> samples[i]"_chr_pos.txt";
                    }
                    prev=$1;
                }
                printf("%d ", $2) >> samples[i]"_chr_pos.txt";
            }
        }
        else {
            print >> samples[i]"_skipped.vcf";
        }
    }
}
END {
    # add newline to file, if it was created
    if (prev != "notset") {
        for (i in samples) {
            printf("\n") >> samples[i]"_chr_pos.txt";
            print samples[i];
        }
    }
}