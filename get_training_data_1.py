import os
import pysam
import math
import numpy as np
from scipy import stats

def get_del_area(vcf_file_path,depth_file_path,chid): #分别获得左断点、右断点的列表
    #输入：vcf文件，深度信息文件，染色体号列表
    open_vcf_file = open(vcf_file_path + chid + '.txt') #vcf每行形式：chr18	10050	10659	1
    vcf_file_list = open_vcf_file.readlines()
    open_depth_file = open(depth_file_path + chid + '_depth.txt')
    depth_file_list = open_depth_file.readlines()
    depth_start = int(depth_file_list[0].strip().split('\t')[1])
    depth_stop = int(depth_file_list[-1].strip().split('\t')[1])
    del_areas = []
    for line in vcf_file_list:
        left_point = int(line.strip().split('\t')[1])
        right_point = int(line.strip().split('\t')[2])
        del_areas.append([max(depth_start,left_point),min(depth_stop,right_point)])
    open_depth_file.close()
    open_vcf_file.close()
    return del_areas , depth_file_list #返回列表，每个元素的type为int,左右扩展50。[[10003, 10709], ……, [78016291, 78016657]]

def get_del_point_inf(vcf_file_path,depth_file_path,chid,window_length=100,stride=75): #获得左或右的范围信息
    #输入：每个位置的深度信息文件，染色体号列表，获得左还是右的信息，粗略筛选bp_range*2范围内是“断点”的断点
    del_areas , depth_file_list = get_del_area(vcf_file_path,depth_file_path,chid) #获得vcf左右断点
    end_point = int(depth_file_list[-1].strip().split('\t')[1])
    final_del_area = []
    for i in range(len(del_areas)): #获得每个左断点和右断点的范围
        del_length = int(del_areas[i][1]) - int(del_areas[i][0]) + 1
        if del_length <= 100:
            final_del_area.append([((int(del_areas[i][1])+int(del_areas[i][0]))//2-49),((int(del_areas[i][1])+int(del_areas[i][0]))//2+50)])
        else:
            for j in range((del_length - window_length)//stride + 1):
                final_del_area.append([int(del_areas[i][0])+j*stride,int(del_areas[i][0])-1+window_length+j*stride])
                if j+1 == ((del_length - window_length)//stride + 1) and ((del_length - (((del_length - window_length)//stride)*stride+window_length)) >= 50) and [int(del_areas[i][1])+1-window_length,int(del_areas[i][1])] != final_del_area[-1]:
                    final_del_area.append([int(del_areas[i][1])+1-window_length,int(del_areas[i][1])])
    print(chid+' get_del_sample_num ... done')
    return final_del_area
