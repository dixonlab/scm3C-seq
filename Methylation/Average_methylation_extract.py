#Description: Calculate average methylation level of each column of file1 on regions from file2
#Both file 1 and 2 have to be modified with Modify_allc_for_merge_for_CpG.py and sorted with modified ID
#Usage: python Average_methylation_extract.py file_with_methylation_level_each_sample_each_column file_target_region output

import sys
ref=open(sys.argv[1],'r')
dfh=open(sys.argv[2],'r')
#rfh=open('Neuron_CTCF_binding_site_all_sample_methylation_+-2kb.txt','w')
rfh=open(sys.argv[3],'w')
data=[]
index=[]
empty=[]
import scipy.stats as stats
rfh.write(dfh.readline().rstrip()+'\t'+'kruskal_pvalue\n')
for a in range(0,9):
        empty.append('NA')
        data.append([])
for q in ref:
        ref_line=q.split()
        ID=ref_line[0]
        target=int(ref_line[0].split('_')[0])
        target_start=int(ref_line[0].split('_')[0])
        target_end=int(ref_line[0].split('_')[1])
        pre=target_start
        if target_start in index:
                a=target_start
                data=map(lambda x:x[index.index(a):],data)
                index=index[index.index(a):]
        if target_start not in index:
                data=[]
                index=[]
                for a in range(0,9):
                        data.append([])
        for i in dfh:
                line=i.split()
                pos=int(line[0].split('_')[0])                               
                if target_start<=pos<=target_end:
                        line=map(lambda x:x.split('(')[0],line)
                        if len(index)!=0:
                                pre=index[-1]
                        if len(index)==0:
                                pre=target_start-1
                        for a in range(pre+1,pos):
                                map(lambda x,y:x.append(y),data,empty)
                                index.append(a)
                        map(lambda x,y:x.append(y),data,line[1:])
                        index.append(pos)
                if pos>target_end:
                        if len(index)!=0:
                                pre=index[-1]
                        if len(index)==0:
                                pre=target_start-1
                        for a in range(pre+1,target_end+1):
                                map(lambda x,y:x.append(y),data,empty)
                                index.append(a)
                        CpG_density_ID=str(len(data[8])-data[8].count('NA'))
                        for q in range(0,4-len(CpG_density_ID)):
                                CpG_density_ID='0'+CpG_density_ID
                        rfh.write(CpG_density_ID+'_'+ref_line[0])
                        temp_datas=[]
                        for a in range(0,14):
                                temp_datas.append([])
                        pp=0
                        for a in data:
                                rfh.write('\t')
                                summ=0
                                coun=0
                                for b in a:
                                        if b!='NA':
                                                coun+=1.0
                                                summ+=float(b)
                                                temp_datas[pp].append(float(b))
                                if coun!=0:
                                        rfh.write(str(summ/coun))
                                if coun==0:
                                        rfh.write('NA')
                                pp+=1
                        size=map(lambda x:len(x),temp_datas)
                        try:
                                kruskal=stats.kruskal(temp_datas[0],temp_datas[1],temp_datas[2],temp_datas[3],temp_datas[4],temp_datas[5],temp_datas[6],temp_datas[7],temp_datas[8],temp_datas[9],temp_datas[10],temp_datas[11],temp_datas[12],temp_datas[13])
                        except:
                                kruskal=['nan','nan']
                        rfh.write('\t'+str(kruskal[1])+'\n')
                        break
