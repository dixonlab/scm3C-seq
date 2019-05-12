import sys                                                                                                                                       
fn=sys.argv[1]
chro=sys.argv[2]
bin_size=sys.argv[3]
ref=open('/pbld/netapp/home/dslee/Refs/mm9/chrNameLength.txt','r')
for i in ref:
        line=i.split()
        if line[0]==chro:
                chro_size=int(line[1])
rfh=open(fn.split('/')[-1]+'_hicrep_matrix_'+chro+'.txt','w')
bin_size=int(bin_size)
dfh=open(fn,'r')
data=[]
for i in range(0,(chro_size/bin_size)+1):
        data.append([])
        for a in range(0,(chro_size/bin_size)+1):
                data[-1].append(0)
for i in dfh:
        line=i.split()
        if line[0].split(':')[0]==chro and line[1].split(':')[0]==chro and 'random' not in i:
                data[int(line[0].split(':')[1].split('-')[0])/bin_size][int(line[1].split(':')[1].split('-')[0])/bin_size]+=int(line[-1])
count=0
for i in range(0,(chro_size/bin_size)+1):
        rfh.write(chro+'\t'+str(i*bin_size)+'\t'+str((i+1)*bin_size))
        for a in range(0,(chro_size/bin_size)+1):
                rfh.write('\t'+str(data[i][a]))
        rfh.write('\n')
