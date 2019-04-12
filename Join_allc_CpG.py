import sys
import gzip
import os
chro_order=['shit']
for i in range(1,24):
        chro_order.append('chr'+str(i))
chro_order.append('chrX')
chro_order.append('chrY')
chro_order.append('chrM')

fn=sys.argv[1]
dfh=gzip.open(fn,'r')
rfh=open(fn.split('/')[-1]+'_CpG_mod_temp.txt','w')
sample=fn.split('_')[4]+'_'+fn.split('_')[9]+'_'+fn.split('_')[10]+'_'+fn.split('_')[11]
rfh.write('0000000000\tSample\t'+sample+'\n')
for i in dfh:
        line=i.split()
        if line[3][:2]=='CG':
                if 'chr'+line[0] in chro_order:
                        if chro_order.index('chr'+line[0])>=10:
                                index=str(chro_order.index('chr'+line[0])*100000000+int(line[1]))
                        if chro_order.index('chr'+line[0])<10:
                                index=str(chro_order.index('chr'+line[0])*1000000000+int(line[1]))
                        rfh.write(index+'\t'+'chr'+'_'.join(line[0:4])+'\t'+str(100*float(line[5])/float(line[6]))+'('+line[5]+','+line[6]+')\n')
rfh.close()
os.system('sort -k 1b,1 '+fn.split('/')[-1]+'_CpG_mod_temp.txt > '+fn.split('/')[-1]+'_CpG_mod.txt')
os.system('rm '+fn.split('/')[-1]+'_CpG_mod_temp.txt')

fn1=sys.argv[1]
fn2=sys.argv[2]


fn=sys.argv[2]
dfh=gzip.open(fn,'r')
dfh=gzip.open(fn,'r')
rfh=open(fn.split('/')[-1]+'_CpG_mod_temp.txt','w')
sample=fn.split('_')[4]+'_'+fn.split('_')[9]+'_'+fn.split('_')[10]+'_'+fn.split('_')[11]
rfh.write('0000000000\tSample\t'+sample+'\n')
for i in dfh:
        line=i.split()
        if line[3][:2]=='CG':
                rfh.write('chr'+'_'.join(line[0:4])+'\t'+str(100*float(line[5])/float(line[6]))+'('+line[5]+','+line[6]+')\n')
rfh.close()
os.system('sort -k 1b,1 '+fn.split('/')[-1]+'_CpG_mod_temp.txt > '+fn.split('/')[-1]+'_CpG_mod.txt')
os.system('rm '+fn.split('/')[-1]+'_CpG_mod_temp.txt')
os.system("join -a1 -a2 -e 'NA' -o 0,1.2,1.3,2.2 "+fn1.split('/')[-1]+'_CpG_mod.txt '+fn2.split('/')[-1]+'_CpG_mod.txt > '+sys.argv[3])
