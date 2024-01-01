import os,sys,subprocess,random
import argparse

## check VCF file for header, otherwise create fake header, this is needed for shapeit2
def add_header(VCF):
	addheader=1;
	File = open(VCF); 
	for line in File: 
		v = line.strip().split('\t'); 
		if line[0] == '#': 
			if v[0] == '#CHROM': addheader=0;
			print >>sys.stdout,line.strip();
		else: 
			if addheader ==1:
				print >>sys.stdout, "#CHROM    POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00000"; 
				print >>sys.stderr, "adding fake header to VCF \n";
				addheader =0;
			alleles = v[4].split(',');
			if len(alleles) > 1: print >>sys.stderr, "filtering variant with multiple alt alleles\n";
			else: print >>sys.stdout,line.strip();
	File.close();


add_header(sys.argv[1]);

