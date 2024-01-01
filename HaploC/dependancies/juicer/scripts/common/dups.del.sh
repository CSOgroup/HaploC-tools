juiceDir=/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/0.Scripts//HaploC/packages/juicer/
site_file=/work/FAC/FBM/DBC/gciriell/default_nonsensitive/Yuanlong/1.Projects/15.Phasing_SV/0.Scripts//HaploC/packages/juicer//restriction_sites/mm10_HindIII.txt
${juiceDir}/scripts/common/fragment.pl SRR7193950_small.PART3_norm.txt test_frag.txt $site_file

awk '$2 == $6 && $4 == $8'  | cut -f1-8 -d" " > intra-fragment.txt


grep chr11 $site_file | awk '{for (i=1; i<=NF; i++) print $i}' > site_chr11.txt


0 chr19 22928985 6400 16 chrY 74451320 21618 60 101M
0 chr19 33321423 9825 16 chr19 33321754 9826 38 101M
0 chr19 36734775 10973 16 chr19 36735018 10974 60 20S81M
0 chr16 13678198 3496 0 chr19 6403094 908 60 101M
16 chr16 10497822 2504 16 chr19 16236101 4064 60 101M
0 chr19 44188206 13216 16 chr19 44188464 13216 60 101M
0 chr19 33014608 9713 0 chrX 108331287 33816 60 13S88M
0 chr19 23645517 6644 16 chr4 25306774 7640 60 101M
16 chr16 71018232 22244 16 chr19 15315740 3746 60 101M
0 chr19 7297383 1073 16 chr19 7297653 1074 60 26S75M


library(data.table)
site_tab = fread('site_chr11.txt')
colnames(site_tab) = 'pos'

frag_tab = fread('intra-fragment.txt')
frag_tab[V2=='chr11']