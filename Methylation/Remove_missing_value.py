#Description: Remove rows with more than X number of 'NA'
#Usage: python Remove_missing_value.py <input> <row with number of NA to remove>
import sys
X=sys.argv[2]
dfh=open(sys.argv[1],'r')
rfh=open(sys.argv[1]+'_All_covered.txt','w')
for i in dfh:
        line=i.split()
        if line.count('NA')>=int(X):
                rfh.write(i)
