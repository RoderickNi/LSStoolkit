# Other_Functions.py


def split_list(lst, n=4):
    '''将一个列表按指定数量合并元素'''
    return [lst[i:i+n] for i in range(0, len(lst), n)]

def count_lines(filename):
    with open(filename, 'r') as file:
        line_count = sum(1 for _ in file)
    return line_count

def GetSeqFromFastq(path):
    '''从fastq文件中提取序列'''
    return [x[1] for x in split_list([x for x in open(path,'r').read().split('\n') if x != ''],n=4)]

complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N': 'N','-':'-'}
def RevComp(dna_sequence):
    '''DNA序列反向互补'''
    return ''.join([complement_dict[base] for base in dna_sequence[::-1]])

def Generate_Subdict(dictionary, n):
    '''利用生成器按批次读取原字典内容，批次数量可自定义'''
    items = list(dictionary.items())
    for i in range(0, len(items), n):
        sub_dict = dict(items[i:i+n])
        yield sub_dict
        
from collections import Counter
def CountListEleNumber(list_A, list_B):
    ''''''
    counter_A = Counter(list_A)
    return [counter_A[x] for x in list_B]

import numpy as np
def Find_Closest(given_array, array_to_compare, num_closest=10):
    '''给定一个形状为（1,2）的数组，在形状为（n,2）的数组中找到欧式距离最近的前n个数组'''
    # 计算欧式距离
    distances = np.sqrt(np.sum((array_to_compare - given_array)**2, axis=1))
    # 获取最小的num_closest个距离的索引
    closest_indexes = np.argpartition(distances, num_closest)[:num_closest]
    # 返回对应的数组元素
    return array_to_compare[closest_indexes]
# 测试
# import matplotlib.pyplot as plt
# given_array = np.array([0.5, 0.5])
# ax.scatter(given_array[0],given_array[1],c='red')
# array_to_compare = np.random.rand(100, 2)
# ax.scatter(array_to_compare[:,0],array_to_compare[:,1],c='black')
# rst=find_closest(given_array, array_to_compare, num_closest=10)
# ax.scatter(rst[:,0],rst[:,1],c='blue',s=100,alpha=0.5)
def Return_Closest_Index(given_array, array_to_compare, num_closest=10):
    '''给定一个形状为（1,2）的数组，在形状为（n,2）的数组中找到欧式距离最近的前n个数组'''
    # 计算欧式距离
    distances = np.sqrt(np.sum((array_to_compare - given_array)**2, axis=1))
    # 获取最小的num_closest个距离的索引
    return np.argpartition(distances, num_closest)[:num_closest]

def remove_common_prefix_suffix(str1, str2):
    i = 0
    while i < min(len(str1), len(str2)) and str1[i] == str2[i]:
        i += 1
    prefix = str1[:i]
    i = 1
    while i <= min(len(str1), len(str2)) and str1[-i] == str2[-i]:
        i += 1
    suffix = str1[-i+1:]
    return str1[len(prefix):-len(suffix)],str2[len(prefix):-len(suffix)]