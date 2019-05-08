#usage: python Modify_supermatrix_for_merge.py input sample_name
import sys                                                                                             
import os
chro_order=['']
for i in range(1,24):
        chro_order.append('chr'+str(i))
chro_order.append('chrX')
chro_order.append('chrY')
chro_order.append('chrM')
fn=sys.argv[1]
dfh=open(fn,'r')
rfh=open(fn.split('/')[-1][:-4]+'_mod.txt','w')
sample1=fn.replace('L2_3','L23').split('_')
sample=sys.argv[2]
rfh.write('000000000000000000000000_Sample\t'+sample+'\n')
for i in dfh:
        line=i.split()
        if line[0] in chro_order and line[3] in chro_order and line[0]==line[3] and line[1]!=line!=[4]:
                if chro_order.index(line[0])>=10:
                        index1=str(chro_order.index(line[0])*10000000000+int(line[1]))
                if chro_order.index(line[0])<10:
                        index1='0'+str(chro_order.index(line[0])*10000000000+int(line[1]))
                if chro_order.index(line[3])>=10:
                        index2=str(chro_order.index(line[3])*10000000000+int(line[4]))
                if chro_order.index(line[3])<10:
                        index2='0'+str(chro_order.index(line[3])*10000000000+int(line[4]))
                index=index1+index2
                rfh.write(index+'_'+'_'.join(line[0:6])+'\t'+line[6]+'\n')
rfh.close()
os.system('sort -k 1b,1 '+fn.split('/')[-1][:-4]+'_mod_temp.txt > '+fn.split('/')[-1][:-4]+'_mod.txt')
os.system('rm '+fn.split('/')[-1]+'_CpG_mod_temp.txt')
