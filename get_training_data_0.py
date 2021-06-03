import pysam
import numpy as np
import os
import math
from scipy import stats

def binarySearch (arr, l, r, x): #二分查找到目标碱基的index
    #depth_file数组，起始，终止，目标数子
    if r >= l: 
        mid = int(l + (r - l)/2)
        if int(arr[mid].strip().split('\t')[1]) == x: 
            return mid 
        elif int(arr[mid].strip().split('\t')[1]) > x: 
            return binarySearch(arr, l, mid-1, x) 
        else: 
            return binarySearch(arr, mid+1, r, x)   
    else: 
        return -1

def merge_del_area(vcf_file_path,chid): #融合有交集的相邻缺失区域，为zero_area函数服务
    open_file_1 = open(vcf_file_path+chid+".txt")
    open_vcf_file = open_file_1.readlines() #'chr18\t10050\t10659\t1\n'
    remp = []
    j = 0
    for i in range(len(open_vcf_file)-1):
        if i == j:
            if (int(open_vcf_file[i+1].strip().split('\t')[1]) - int(open_vcf_file[i].strip().split('\t')[2])) > 0:
                remp.append(int(open_vcf_file[i].strip().split('\t')[1]))
                remp.append(int(open_vcf_file[i].strip().split('\t')[2]))
                j += 1
            elif (int(open_vcf_file[i+1].strip().split('\t')[1]) - int(open_vcf_file[i].strip().split('\t')[2])) <= 0:
                remp.append(int(open_vcf_file[i].strip().split('\t')[1]))
                if open_vcf_file[i+1].strip().split('\t')[2] >= open_vcf_file[i].strip().split('\t')[2]:
                    remp.append(int(open_vcf_file[i+1].strip().split('\t')[2]))
                else:
                    remp.append(int(open_vcf_file[i].strip().split('\t')[2]))
                j += 2
    merged_del_area = np.array(remp).reshape((-1,2))
    print(chid+' merge_del_area ... done')
    open_file_1.close()
    return merged_del_area #numpy(-1,2)


def zero_area(vcf_file_path,depth_file_path,chid): #得到最终采用的未发生变异的区域，为get_zero_feature服务
    zero_area = []
    del_zero_area_index = []
    merged_del_area = merge_del_area(vcf_file_path,chid) #已经融合好的变异区域
    merged_del_area_len = len(merged_del_area) #融合好的变异区域个数
    open_depth = open(depth_file_path+chid+"_depth.txt") #'chr18\t10003\t1\n'
    open_depth_file = open_depth.readlines()
    open_depth_file_len = len(open_depth_file) #记录深度文件的每行信息和len
    for i in range(merged_del_area_len):
        l_del_point = merged_del_area[i][0]
        r_del_point = merged_del_area[i][1] #记录融合好的左右断点
        l_point_index = binarySearch(open_depth_file,0,open_depth_file_len - 1,l_del_point)
        r_point_index = binarySearch(open_depth_file,0,open_depth_file_len - 1,r_del_point) #找到融合好的左右断点在深度文件中的索引                                  
        if l_point_index != -1 and r_point_index != -1 and (l_point_index-r_point_index) == (l_del_point-r_del_point):
            zero_area.append(l_del_point)
            zero_area.append(l_point_index)
            zero_area.append(r_del_point)
            zero_area.append(r_point_index)
    #if (zreo_area[0]-int(open_depth_file[0].strip().split('\t')[1])) == 
#     if (zero_area[0]-int(open_depth_file[0].strip().split('\t')[1])) == (zero_area[1]-binarySearch(open_depth_file,0,open_depth_file_len - 1,int(open_depth_file[0].strip().split('\t')[1]))):
#         zero_area.insert(0,binarySearch(open_depth_file,0,open_depth_file_len - 1,int(open_depth_file[0].strip().split('\t')[1])))
#         zero_area.insert(0,int(open_depth_file[0].strip().split('\t')[1]))
#     if (zero_area[-2]-int(open_depth_file[-1].strip().split('\t')[1])) == (zero_area[-1]-binarySearch(open_depth_file,0,open_depth_file_len - 1,int(open_depth_file[-1].strip().split('\t')[1]))):
#         zero_area.append(int(open_depth_file[-1].strip().split('\t')[1]))
#         zero_area.append(binarySearch(open_depth_file,0,open_depth_file_len - 1,int(open_depth_file[-1].strip().split('\t')[1])))
    zero_area.insert(0,binarySearch(open_depth_file,0,open_depth_file_len - 1,int(open_depth_file[0].strip().split('\t')[1])))
    zero_area.insert(0,int(open_depth_file[0].strip().split('\t')[1]))
    zero_area.append(int(open_depth_file[-1].strip().split('\t')[1]))
    zero_area.append(binarySearch(open_depth_file,0,open_depth_file_len - 1,int(open_depth_file[-1].strip().split('\t')[1])))
    #print(zero_area[:5],len(zero_area))
    zero_array = np.array(zero_area).reshape((-1,4))[:,[0,2]]
    for i in range(len(zero_array)):
        if (zero_array[i][1]-zero_array[i][0])<100000: #10w bp哪没有发生缺失，记为可用范围
            del_zero_area_index.append(i)
    f_zero_array = np.delete(zero_array,del_zero_area_index,axis=0)
    final_zero_array = []
    for i in range(len(f_zero_array)):
        l_zero_point = f_zero_array[i][0]
        r_zero_point = f_zero_array[i][1]
        l_zero_point_index = binarySearch(open_depth_file,0,open_depth_file_len - 1,l_zero_point)
        r_zero_point_index = binarySearch(open_depth_file,0,open_depth_file_len - 1,r_zero_point)
        if (l_zero_point_index-r_zero_point_index) == (l_zero_point-r_zero_point):
            final_zero_array.append(f_zero_array[i])
    open_depth.close()
    print(chid+' zero_area ... done')
    return np.array(final_zero_array),open_depth_file #numpy(-1,2)


def get_zero_feature_num(vcf_file_path,depth_file_path,chid,get_num=10000): #得到未发生变异区域的特征
    #zero_area返回值。get_num在给定染色体中得到给定个数的未发生变异区域的特征
    #final_zero_array = array
    final_zero_array,open_depth_file = zero_area(vcf_file_path,depth_file_path,chid)
    final_zero_array_length = final_zero_array[:,1] - final_zero_array[:,0]
    final_zero_array_length_sum = np.sum(final_zero_array,axis=0)[1] - np.sum(final_zero_array,axis=0)[0]
    every_zero_area_num = (final_zero_array_length / final_zero_array_length_sum * get_num)+1
    print(chid+' get_zero_sample_num ... done')
    return final_zero_array,every_zero_area_num.astype(int),open_depth_file,len(open_depth_file)

#final_zero_array:array[[  109692,   393327],
#                        [  393436,   512241],
#                       [  667520,   837074],
#                        ...
#                        [77831572, 77964969]])

#every_zero_area_num:array([ 39,  16,  23,...18])



def get_final_zero_area(final_zero_array,every_zero_area_num,reserved_length=1000,window_length=100): #返回生成为变异区域的起始位点
    random = np.random.RandomState(4396)  #随机数种子，相同种子下每次运行生成的随机数相同
    final_zero_areas = []
    for i in range(len(final_zero_array)):
        single_area_zero_num = random.randint(low=final_zero_array[i][0]+reserved_length,high=final_zero_array[i][1]-reserved_length, size=every_zero_area_num[i])
        for j in range(len(single_area_zero_num)):
            final_zero_areas.append([single_area_zero_num[j],single_area_zero_num[j]-1+window_length])
    return final_zero_areas#!/usr/bin/env python

