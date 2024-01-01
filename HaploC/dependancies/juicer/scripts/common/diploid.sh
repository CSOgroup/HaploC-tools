#!/bin/bash

# stem=/aidenlab/work/neva/patski_


merged_nodups_file=/cso/gscratch/Yuanlong/1.Data/public_hic_data/T47D/merge_bp/aligned/merged_nodups.txt  
chr_pos_file=unknown_chr_pos.txt
SNP_file=unknown_paternal_maternal.txt

juiceDir=/cso/gred/Yuanlong/3.Software_packages/juicer/scripts/common/

# merged_nodups_file=/mnt/ndata/Yuanlong/1.Projects/15.Phasing_SV/1.Data/GM12878/GSM1551550_HIC001_merged_nodups.txt

for i in chr18 chr19
  do
  awk -v chr=$i '$2==chr && $6==chr && $9 >= 10 && $12 >= 10' $merged_nodups_file | ${juiceDir}/diploid.pl -s $chr_pos_file -o $SNP_file > diploid_${i}.txt
done



cat diploid_chr1.txt diploid_chr10.txt diploid_chr11.txt diploid_chr12.txt diploid_chr13.txt diploid_chr14.txt diploid_chr15.txt diploid_chr16.txt diploid_chr17.txt diploid_chr18.txt diploid_chr19.txt diploid_chr2.txt diploid_chr21.txt diploid_chr22.txt diploid_chr3.txt diploid_chr4.txt diploid_chr5.txt diploid_chr6.txt diploid_chr7.txt diploid_chr8.txt diploid_chr9.txt diploid_chrX.txt diploid_chrY.txt > diploid.txt

awk -f ${juiceDir}/scripts/diploid_split.awk diploid.txt 
sort -k2,2d -m maternal.txt maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > mat.txt 
${juiceDir}/scripts/juicer_tools pre mat.txt maternal.hic hg19 
sort -k2,2d -m paternal.txt paternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > pat.txt 
/aidenlab/juicebox pre pat.txt paternal.hic hg19

