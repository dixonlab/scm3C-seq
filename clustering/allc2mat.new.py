import sys
import gzip
import numpy as np
from multiprocessing import Pool

def add(x,mode):
	if mode==1:
		if x[5]>2 or x[3][1]=='N':
			y=[0,0,0,0]
		elif x[3][1]!='G':
			y=[x[4],x[5],0,0]
		else:
			y=[0,0,x[4],x[5]]
		return y
	if mode==0:
		y=[0 for i in range(6)]
		if x[5]>2:
			return y
		if x[3]=='CCC':
			y[0]+=x[4]
			y[1]+=x[5]
		if x[3][1]=='G':
			y[4]+=x[4]
			y[5]+=x[5]
		elif x[3][1]!='N':
			y[2]+=x[4]
			y[3]+=x[5]
		return y

def mapping(chr):
	x = bin[bin[:,0]==('chr'+chr), 1:3].astype(int)
	if len(x)==0:
		return
	fin = gzip.open(sys.argv[1], 'r')
	fout = open(out_folder + chr + '.' + sys.argv[3] + '.mC.txt', 'w')
	ans0 = [0 for i in range(6)]
	if 'chr'+chr in idx.keys():
		chr = 'chr' + chr
	elif not (chr in idx.keys()):
		fin.close()
		fout.close()
		return(ans0)
	fin.seek(idx[chr])
	print(chr)
	tmp = [chr]+[-1 for i in range(6)]
	y = [tmp[:]]
	for i in range(len(x)):
		ans = [0,0,0,0]
		if len(y)>0 and y[-1][1]>x[i,1]:
			j = len(y)
			while j-1>=0 and y[j-1][1]>x[i,1]:
				j -= 1
			while j-1>=0 and y[j-1][1]>x[i,0]:
				j -= 1
				result = add(y[j],1)
				ans = [ans[k] + result[k] for k in range(4)]
			y = y[j:]
		elif len(y)>0 and y[-1][1]>x[i,0]:
			j = len(y)
			while j-1>=0 and y[j-1][1]>x[i,0]:
				j -= 1
				result = add(y[j],1)
				ans = [ans[k] + result[k] for k in range(4)]
			y = y[j:]
		elif len(tmp)==7 and tmp[0]==chr:
			y = []
			while tmp[1]<=x[i,0]:
				tmp = fin.readline().decode('utf-8').strip().split('\t')
				if len(tmp)<7 or tmp[0]!=chr:
					break
				tmp[1], tmp[4], tmp[5] = int(tmp[1]), int(tmp[4]), int(tmp[5])
				result = add(tmp, 0)
				ans0 = [ans0[k] + result[k] for k in range(6)]
			if len(tmp)==7 and tmp[0]==chr:
				if tmp[1]<=x[i,1]:
					result = add(tmp, 1)
					ans = [ans[k] + result[k] for k in range(4)]
				y.append(tmp)
		if len(tmp)==7 and tmp[0]==chr:
			while tmp[1]<=x[i,1]:
				tmp = fin.readline().decode('utf-8').strip().split('\t')
				if len(tmp)<7 or tmp[0]!=chr:
					break
				tmp[1], tmp[4], tmp[5] = int(tmp[1]), int(tmp[4]), int(tmp[5])
				if tmp[1]<=x[i,1]:
					result = add(tmp,1)
					ans = [ans[k] + result[k] for k in range(4)]
				result = add(tmp, 0)
				ans0 = [ans0[k] + result[k] for k in range(6)]
				y.append(tmp)
		if ans[1]>0 and ans[3]>0:
			result = [float(ans[0])/ans[1], float(ans[2])/ans[3]]
		elif ans[1]>0 and ans[3]==0:
			result = [float(ans[0])/ans[1], 0.0]
		elif ans[1]==0 and ans[3]>0:
			result = [0.0, float(ans[2])/ans[3]]
		else:
			result = [0.0, 0.0]
		fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(ans[0],ans[1],result[0],ans[2],ans[3],result[1]))
	while sys.argv[3]!='test' and len(tmp)==7 and tmp[0]==chr:
		tmp = fin.readline().decode('utf-8').strip().split('\t')
		if len(tmp)<7 or tmp[0]!=chr:
			break
		tmp[1], tmp[4], tmp[5] = int(tmp[1]), int(tmp[4]), int(tmp[5])
		result = add(tmp,0)
		ans0 = [ans0[k] + result[k] for k in range(6)]
	fin.close()
	fout.close()
	return(ans0)

sample = sys.argv[2].split('/')[-1]
ncpus = int(sys.argv[4])
global bin,out_folder,idx
idx = {}
out_folder = sys.argv[2] + '_'
fin = open(sys.argv[1]+'.idx')
for line in fin:
	tmp = line.strip().split()
	if len(tmp)!=2:
		break
	idx[tmp[0]] = int(tmp[1])
fin.close()
chr = ['X']+[str(i+1) for i in range(22)]
if sys.argv[3][-3:]=='bin':
	bin = np.load('/gale/netapp/home/zhoujt/genome/hg19/hg19.' + sys.argv[3] + '.npy')
elif sys.argv[3]=='gene':
	bin = np.load('/gale/netapp/home/zhoujt/resource/gene/hg19/GENCODE/gencode.v29.npy')
elif sys.argv[3]=='test':
	bin = np.loadtxt(sys.argv[5],dtype=np.str)
if bin.dtype.str[1]=='S':
	bin = np.core.defchararray.decode(bin)
p = Pool(ncpus)
result = p.map(mapping,chr)
p.close()
if sys.argv[3]!='test':
	result = np.sum(np.array(result),axis=0)
	fout = open(sys.argv[2] + '.' + sys.argv[3] + '.tot.txt', 'w')
	fout.write(sample + '\t')
	fout.write('\t'.join([str(result[i]) for i in range(6)])+'\t')
	fout.write('{0}\t{1}\t{2}\n'.format(float(result[0])/result[1],float(result[2])/result[3],float(result[4])/result[5]))
	fout.close()

