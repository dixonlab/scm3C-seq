import sys
import os


fn1=sys.argv[1].split('/')[-1]
fn2=sys.argv[2].split('/')[-1]
dfh1=open(fn1,'r')
dfh2=open(fn2,'r')
line1=dfh1.readline().split()
line2=dfh2.readline().split()
o=[]
count=0
for i in line1[1:]:
        o.append(str(count+2))
        count+=1
command="join -a1 -a2 -e 'NA' -o 0"
for i in o:
        command+=',1.'+i
for i in o:
        command+=',2.'+i
print command
dfh1.close()
dfh2.close()
os.system(command+' '+fn1+' '+fn2+' > '+sys.argv[3])
os.system('rm '+fn1+' '+fn2)
