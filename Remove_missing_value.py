import sys
dfh=open(sys.argv[1],'r')
rfh=open(sys.argv[1]+'_All_covered.txt','w')
for i in dfh:
        line=i.split()
        if line.count('NA')==0:             
                rfh.write(i)
