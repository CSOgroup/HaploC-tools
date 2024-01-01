import os,sys,subprocess,random,math
from mstkruskal import Graph

def pairwiseLD(haplist,nsamples,i,j):
	counts = [0,0];
	for r in range(nsamples): 
		if haplist[r][i] == '0' and haplist[r][j] == '0': counts[0] +=1; 
		elif haplist[r][i] == '1' and haplist[r][j] == '1': counts[0] +=1; 
		elif haplist[r][i] == '0' and haplist[r][j] == '1': counts[1] +=1; 
		elif haplist[r][i] == '1' and haplist[r][j] == '0': counts[1] +=1; 
	return counts;

## read the file with list of sampled haplotypes and print summary pairwise statistics
def summarize_haps(haplist,varindexes,THRESH=0.9):
	nsamples = len(haplist); 
	v = len(haplist[0]); 
	edges =[];
	W = 5;
	size = v;

	discarded=0;
	for i in range(size):
		for j in range(i+1,min(v,i+W)):
			counts = pairwiseLD(haplist,nsamples,i,j);
			f = float(counts[0]+0.5)/(counts[0]+counts[1]+1);
			if f > THRESH: 
				edges.append([i,j,f,'00']);  
				#print edges[-1]
			elif 1.0-f > THRESH: 
				edges.append([i,j,1.0-f,'01']); 
				#print edges[-1]
		#print i,j,counts;

	g= Graph(size+W);
	g.addEdgeList(edges);
	mst_edges = g.KruskalMST();
	mst_edges = sorted(mst_edges,key=lambda item: item[0])
	fid = 1;
	for edge in mst_edges: 
		blocks =2;
		## (1-q)^2 + q^2 = edge[3],  (1-edge[3])/2 = q(1-q) 
		## q =   0.5 + sqrt (1-4t)/2
		## for edge[3]= 0.99, q = 0.9947
		t = (1.0-float(edge[2]))/2; 
		q = 0.5 + math.sqrt(1-4*t)/2; 
		qv = chr(int(-10*math.log(1-q,10))+33);
		print >>sys.stdout,"2 FRAGMENT_" + str(fid),'0 -1 -1',varindexes[edge[0]],edge[3][0],varindexes[edge[1]],edge[3][1],qv+qv;
		#print edge[0],edge[1],edge[2],edge[3],varindexes[edge[0]],varindexes[edge[1]];
		#print 'stats',edge[2],t,q,qv
		fid +=1;
	print >>sys.stderr,"edges",len(mst_edges),size;


def read_shapeit_samples(HAPfile):
	File = open(HAPfile,'r'); 
	nsamples=0; haplist = []
	for line in File: 
		
		if nsamples ==0: ## first line in new format has variant indexes
			varindexes = map(int, line.strip().split());  
		else: 
			haplist.append(line.strip()); 
		nsamples +=1; 
	File.close(); 
	nsamples -=1; 
	print >>sys.stderr, "read",nsamples,'haplotypes from file'; 
	return [haplist,varindexes]

########################################################################################################

HAPfile = sys.argv[1]
[haplist,varindex] = read_shapeit_samples(HAPfile);
if len(sys.argv) > 2: summarize_haps(haplist,varindex,float(sys.argv[2]));
else: summarize_haps(haplist,varindex);



