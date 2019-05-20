import sys
from collections import Counter
import numpy as np

infile = sys.argv[1]
res = int(sys.argv[3])
fin = open(infile, 'r')
if res>=1000000:
	k = str(res/1000000)+'mb'
elif res>=1000:
	k = str(res/1000)+'kb'

if sys.argv[4][:2]=='hg':
	chrom = [str(i+1) for i in range(22)] + ['X']
else:
	chrom = [str(i+1) for i in range(19)] + ['X']

count = Counter()
for line in fin:
	tmp = line.strip().split('\t')
	tmp[1], tmp[3] = tmp[1][3:], tmp[3][3:]
	if (np.abs(int(tmp[4])-int(tmp[2])) >= 10000) and (tmp[1]==tmp[3]) and (tmp[1] in chrom):
		tmp[2], tmp[4] = int(tmp[2]) // res, int(tmp[4]) // res
		if tmp[2]>tmp[4]:
			tmp[2], tmp[4] = tmp[4], tmp[2]
		if tmp[2]!=tmp[4]:
			pos = '-'.join([tmp[1], str(tmp[2]), str(tmp[4])])
			count[pos] += 1

fin.close()
fout = {c:open(sys.argv[2] + '_chr' + c + '.txt', 'w') for c in chrom}
for key in count:
	tmp = key.split('-')
	fout[tmp[0]].write('{0}\t{1}\t{2}.0\n'.format(tmp[1], tmp[2], count[key]))

for c in chrom:
	fout[c].close()


