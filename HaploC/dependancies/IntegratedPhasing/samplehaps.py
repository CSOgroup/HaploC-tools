import os,sys,subprocess,random
import argparse

#################################

SHAPEIT = "%s/%s" % ( sys.argv[4],  "/HaploC/dependancies/SHAPEIT2/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit")

## first line of output file has list of variant-ids (relative to VCF file) 
def sample_haplotypes(VCF,out_prefix,nsamples):
	HAPfile = out_prefix + '.hapsamples'; 
	snptable = {}; variantid=0;

	VCF_file = open(VCF);
	for line in VCF_file:
		if line[0]== '#': continue;
		var = line.split();
		chrom = var[0]; position = int(var[1]); allele1 = var[3]; alleles = var[4].split(','); allele2 = alleles[0];
		genotype = var[9].split(':');
		variantid +=1;
		snptable[(chrom,position,allele1)] = [genotype[0],1,variantid]; 
	VCF_file.close();


	if nsamples > 1: F1 = open(HAPfile,'w'); # new output file

	for r in range(nsamples):
		# seed = random.randint(1,10000000);
		seed = r;
		outfile = out_prefix + '.phased'; logfile = out_prefix + '.shapeit.log';
		# samplehap = SHAPEIT + ' -convert' + ' --seed ' + `r` + ' --input-graph ' + VCF + '.graph' + ' --output-sample ' + outfile + ' -L ' + logfile;
		samplehap = SHAPEIT + ' -convert' + ' --seed ' + `r` + ' --input-graph ' + out_prefix + '.graph' + ' --output-sample ' + outfile + ' -L ' + logfile + ' > /dev/null';		
		subprocess.call(samplehap,shell=True);
		hap1 = []; hap2 = [];
		File = open(outfile + '.haps','r');
		for line in File: 
			var = line.strip().split(); h1 = var[5]; h2 = var[6]; 
			if h1 == h2: continue; ## homozoygous variant
			if r ==0: 
				varid = snptable[(var[0],int(var[2]),var[3])][2]
				print >>F1, varid, ## variant-id 
			hap1.append(h1); hap2.append(h2); 
		File.close();
		if r ==0: print >>F1,'\n',
		if nsamples > 1: print >>F1, ''.join(hap1); ## only write one of the two haplotypes
		if (r+1)%500==0: print >>sys.stderr,'iter',r+1;
	if nsamples > 1: F1.close(); 
	subprocess.call('rm -f ' + outfile + '.* ' + logfile,shell=True);


if len(sys.argv) < 3: 
	print >>sys.stderr, "program to sample haplotypes from SHAPEIT2 HMM ";
	print >>sys.stderr, "set path to shapeit before running this program, the shapeit graph file also needs to exist ";
	print >>sys.stderr, "usage: python samplehaps.py VCFfile number_haps(integer)";
	sys.exit();
sample_haplotypes(sys.argv[1],sys.argv[2],int(sys.argv[3])); 

