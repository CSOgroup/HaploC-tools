

#### Requirements: 

1. HapCUT2: https://github.com/vibansal/HapCUT2
2. SHAPEIT2: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html 
3. 1000 Genomes reference panel haplotypes: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
4. BAM file and VCF file for individual to be phased 

#### Notes: 

The code was developed and tested using python2.7. Currently, it doesn't work for python versions 3 or more. 

'VCF' refers to the path to the input VCF file

#### Steps: 

1. Run extractHAIRS on the BAM and VCF file  to obtain Hi-C fragments (see HapCUT2 instructions for how to do this, use option --hic 1)

2. Run SHAPEIT2 on the VCF file using the 1000 Genomes reference panel:

```
shapeit -check --input-vcf VCF -R reference_panel.hap.gz reference_panel.legend.gz reference_panel.samples --output-log VCF.log
shapeit --input-vcf VCF -R reference_panel.hap.gz reference_panel.legend.gz reference_panel.samples -M genetic_map.txt --output-graph VCF.graph --exclude-snp VCF.check.snp.strand.exclude
```
The files reference_panel.hap.gz, reference_panel.legend.gz and reference_panel.samples (for each chromosome) can be downloaded from the link above. 

3. use samplehaps.py to sample 'N' haplotype pairs for the individual (we have used N=1000). The program will output the samples to the file VCF.hapsamples where VCF is the path to the input VCF file

```
python samplehaps.py VCF 1000 
```

4. use encodereads.py to represent the haplotype samples as "pseudo-reads" 

```
python encodereads.py VCF.hapsamples > sample.pseudo_reads
```

5. concatenate the "Hi-C fragments" file and the sample.pseudo_reads to obtain the new "fragment file"

6. run HapCUT2 on the new "fragment file" to obtain phased haplotypes 
