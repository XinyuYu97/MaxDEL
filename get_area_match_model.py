#此文件是为了得到该area内的，满足阈值长度的del个数，软切个数。是从画图的程序中截取修改来的
import sys
import pysam
import numpy as np
import os
from PIL import Image
import math

#超级憨批的pysam库记录
##query_sequence方法：返回当前read的序列，此序列包括MIS字段，因为D字段是删除，没有碱基信息
##reference_start和reference_end方法：返回的是read匹配上的ref位点。因为I字段匹配不上，设置为折叠I字段；S字段为软切，在别的地方有匹配
# 所以reference_end-reference_start=M+D(可能会有1的差距)

def read_cigartuples_remove_redundancy(read_cigartuples):
    read_cigartuples_1 = [read_cigartuples[0]]
    for i in range(1,len(read_cigartuples)):
        if read_cigartuples[i][0] == read_cigartuples_1[-1][0]:
            read_cigartuples_1[-1][1] += read_cigartuples[i][1]
        else:
            read_cigartuples_1.append(read_cigartuples[i])
    return np.array(read_cigartuples_1)

#def get_area_read_name(l_point,r_point,bam_file_path):
def get_area_read_seq_inf(open_bam_file,chid,start,stop):
    #输入的bam_file_path为单独染色体，使用该方法时需要将bam文件按照染色体切分
    #M I D S，即：0 1 2 4
    read_temp,cigartuples_temp = [],[]
    for read in open_bam_file.fetch(contig=chid[3:],start=max(0,start-10000),stop=stop+10000):
        #此处的start和stop应为cvf记录的断点。并且包括范围应该比vcf所记录的长度更大一点
        
        if read.is_secondary == False: #筛选只有MIDS的read（几乎read也只有这几个）
            read_cigartuples = np.array(read.cigartuples) #读取cigar字段为:[n*(匹配类型，匹配长度)]
            maped_type = read_cigartuples[:,0]
            if (3 not in maped_type[1:-1]) and (4 not in maped_type[1:-1]) and (5 not in maped_type[1:-1]) and (6 not in maped_type[1:-1]) and (7 not in maped_type[1:-1]) and (8 not in maped_type[1:-1]):
                if (read_cigartuples[0][0] in [0,2,4]) and (read_cigartuples[-1][0] in [0,2,4]):
                    read_name = read.query_name
                    if int(read_cigartuples[0][0]) != 4:
                        read_start = read.reference_start + 1
                    else:
                        read_start = read.reference_start - read_cigartuples[0][1] + 1
                    if int(read_cigartuples[-1][0]) != 4:
                        read_end = read.reference_end
                    else:
                        read_end = read.reference_end + read_cigartuples[-1][1]
                    #以上计算完可用read的起始和终止
            if (read_end<start) or (read_start>stop):
                continue
            else:
                read_temp.append([read_name,read_start,read_end])
#-----------处理掉长度小于3匹配部分,变为为M，变更长于3的insert为1个像素，不需要就直接删掉---------------------------------
                
                #在cigartuples中，M=0，I=1，D=2，N=3，S=4，H=5，P=6，使用M I D S，即：0 1 2 4
                del_mini_area_index = []
                for j in range(len(read_cigartuples)): #删除小于等于3的 I
                    if (read_cigartuples[j][1] <= 3) and (read_cigartuples[j][0] == 1):
                        del_mini_area_index.append(j)
                read_cigartuples = np.delete(read_cigartuples,del_mini_area_index,axis=0)
                for k in range(len(read_cigartuples)): #剩余小于等于3的都变为 M
                    if read_cigartuples[k][1] <= 3:
                        read_cigartuples[k][0] = 0
                        
                read_cigartuples = read_cigartuples_remove_redundancy(read_cigartuples)
                
                for i in range(1,len(read_cigartuples)-1):
                    if read_cigartuples[i][0] == 1:
                        if read_cigartuples[i-1][1] >= read_cigartuples[i+1][1]:
                            read_cigartuples[i-1][1] -= 1
                        else:
                            read_cigartuples[i+1][1] -= 1
                        read_cigartuples[i][1] = 1
                del_model_index = []
                for u in range(len(read_cigartuples)):
                    if read_cigartuples[u][1] == 0:
                        del_model_index.append(u)
                read_cigartuples = np.delete(read_cigartuples,del_model_index,axis=0)
                read_cigartuples = read_cigartuples_remove_redundancy(read_cigartuples)
#                 for l in range(1,len(read_cigartuples)-1):
#                     if (read_cigartuples[l][1] <= 3) and (read_cigartuples[l][0] != 1):
#                         if read_cigartuples[l-1][1] >= read_cigartuples[l+1][1]:
#                             read_cigartuples[l][0] = read_cigartuples[l-1][0]
#                         else:
#                             read_cigartuples[l][0] = read_cigartuples[l+1][0]
#-----------处理掉长度小于3匹配部分,变为为M，变更长于3的insert为1个像素，不需要就直接删掉---------------------------------
                cigartuples_temp.append(read_cigartuples)
    return read_temp,cigartuples_temp


def start_cigar_stop(read_temp,cigartuples_temp,start,stop):
#得到cigar在单独小区域内的读数
    zero_threshold = []
    for i in range(len(read_temp)):
        temp = []
        read_start = read_temp[i][1]
        read_stop = read_temp[i][2]
        #read起始小于等于候选区域
        if read_start >= start:
            for j3 in range(len(cigartuples_temp[i])):
                temp.append([cigartuples_temp[i][j3][0],(min((read_start+cigartuples_temp[i][j3][1]-start),(stop-start+1))-(read_start-start))])
                read_start += cigartuples_temp[i][j3][1]
                if read_start > stop:
                    break
            zero_threshold.append(temp)
            
        elif read_start < start:
            #得到起始的read index
            for k2 in range(len(cigartuples_temp[i])):
                if (read_start + cigartuples_temp[i][k2][1]) >= start:
                    read_start_index = read_start
                    cigartuples_index = k2
                    break
                else:
                    read_start += cigartuples_temp[i][k2][1]
            for k4 in range(cigartuples_index,len(cigartuples_temp[i])):
                temp.append([cigartuples_temp[i][k4][0],(min(read_start_index+cigartuples_temp[i][k4][1]-start,stop-start+1)-max(0,read_start_index-start))])
                read_start_index += cigartuples_temp[i][k4][1]
                if read_start_index > stop:
                    break
            zero_threshold.append(temp)
    return zero_threshold,len(zero_threshold)

def count_match_num(area_match_list,read_num):
    count_del,count_soft_clip,count_primary = 0,0,0
    for i in range(len(area_match_list)):
        for j in range(len(area_match_list[i])):
            if area_match_list[i][j][0] == 2:
                if area_match_list[i][j][1] >= 15:
                    count_del += (area_match_list[i][j][1]//15)
            elif area_match_list[i][j][0] == 4:
                count_soft_clip += 1
                
        for k in range(len(area_match_list[i])): #count count_primary
            if area_match_list[i][k][0] in [0,1,2,3,5]:
                count_primary += 1
                break
    return count_del,count_soft_clip,count_primary,read_num

def match_model(open_bam_file,chid,start,stop):
    #open_bam_file = pysam.AlignmentFile(bam_file_path+'sorted_final_merged_'+chid+'.bam','rb')
    read_temp,cigartuples_temp = get_area_read_seq_inf(open_bam_file,chid,start,stop)
    zero_threshold,read_num = start_cigar_stop(read_temp,cigartuples_temp,start,stop)
    count_del,count_soft_clip,count_primary,read_num = count_match_num(zero_threshold,read_num)
    return count_del,count_soft_clip,count_primary,read_num