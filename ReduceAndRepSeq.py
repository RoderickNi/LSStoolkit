# ReduceAndRepSeq.py
# 输入：fasta格式文件
# 输出：降维去冗余后获得的代表序列

import sys
import time
import umap
import pickle
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction.text import CountVectorizer
from IO import Read_Fasta
from ForSklearn import Get_DBSCAN_CenterPosition



def Return_Closest_Index(given_array, array_to_compare, num_closest=10):
    '''给定一个形状为（1,2）的数组，在形状为（n,2）的数组中找到欧式距离最近的前n个数组'''
    # 计算欧式距离
    distances = np.sqrt(np.sum((array_to_compare - given_array)**2, axis=1))
    # 获取最小的num_closest个距离的索引
    return np.argpartition(distances, num_closest)[:num_closest]

def KmerFrequencyStatistics(Fasta_Path,Min_ngram=1,Max_ngram=6):
    '''
    Vectorization of DNA sequence features by kmer method
    Parameter Detailed Explanation
    '''
    print(f'Kmer words frequency statistics are being performed on all sequences.')
    print(f'The range of n-gram is from {str(Min_ngram)} to {str(Max_ngram)}.')
    NameLst,SeqLst=Read_Fasta(Fasta_Path).Get_NameLst_SeqLst()   # 从Fasta文件中读取 DNA 序列数据
    vectorizer = CountVectorizer(analyzer='char',ngram_range=(Min_ngram, Max_ngram)) # 创建 CountVectorizer 对象，并设置 ngram_range 参数为 (1, 1)
    arr = vectorizer.fit_transform(SeqLst).toarray() # 对 DNA 序列进行向量化处理
    print('Run successfully!')
    return NameLst,arr

def DimensionalReductionVectors(SeqArr,preprocess='Standard',method='PCA_then_UMAP',n_components=2,random_state=0,min_dist=0.0):
    ''''''
    print(f'Got {SeqArr.shape[0]} sequences.\nThe number of features is {SeqArr.shape[1]}.')
    print(f'Data preprocessing is currently underway, and the preprocessing method is {preprocess}.')
    if  preprocess == 'Standard':
        SeqArr = StandardScaler().fit_transform(SeqArr)
    elif  preprocess =='Robust':
        SeqArr = RobustScaler().fit_transform(SeqArr)    
    print(f'The dimensionality reduction is being performed on the preprocessed data, and the dimensionality reduction algorithm is {method}.')
    if method == 'PCA':
        reduced_arr = PCA(n_components=n_components,
                           random_state=random_state).fit_transform(SeqArr)    
    elif method == 'UMAP':
        reduced_arr = umap.UMAP(n_neighbors=50, 
                                 min_dist=min_dist,
                                 metric='cosine',
                                 n_components=n_components).fit_transform(SeqArr)
    elif method == 'PCA_then_UMAP':
        reduced_arr = umap.UMAP(n_neighbors=50, 
                                 min_dist=min_dist,
                                 metric='cosine',
                                 n_components=n_components).fit_transform(PCA(n_components=100,
                                                                              random_state=random_state).fit_transform(SeqArr))
    print('Run successfully!')
    return reduced_arr

def GetRepSeq(NameLst,SeqArr,SampleSize,Fasta_Path,Otpt_RepSeq,radius=1,num_closest=5):
    print('Determine the cluster of clusters based on density clustering, and calculate the centroid coordinates of each cluster.')
    SeqArr=MinMaxScaler(feature_range=(-50,50)).fit_transform(SeqArr)# 数据预处理缩放至-50 ~ 50 （100 X 100）
    min_samples,eps=int(SeqArr.shape[0]*(1/(SampleSize*2))),radius # 每个簇最少支持的点的个数，和同簇点识别半径
    dbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(SeqArr) # 创建DBSCAN模型并拟合
    centers=Get_DBSCAN_CenterPosition(SeqArr,dbscan)
    print('Extract representative sequences of each cluster based on the centroid coordinates of each cluster.')
    NameSeq_dict=Read_Fasta(Fasta_Path).GetNameSeq()
    with open(Otpt_RepSeq,'w',encoding='utf-8') as f:
        pass
    for cluster,INFO in centers.items():
        position,total_num=INFO['POS'],INFO['NUM']
        closest_indexes = Return_Closest_Index(position,SeqArr,num_closest=num_closest)
        for seqname in [NameLst[x] for x in closest_indexes]:
            print(cluster,seqname,position,total_num)
            with open(Otpt_RepSeq,'a',encoding='utf-8') as f:
                print(f'>Allele-{cluster}@{seqname}::{str(total_num)}\n{NameSeq_dict[seqname]}',file=f)
    print('Run successfully!')
    
    
if __name__ == '__main__':
    
    # 参数默认值
    Min_ngram,Max_ngram = 1,5
    SampleSize=50
    preprocess,method='Standard','PCA_then_UMAP'
    n_components=2
    random_state=0
    min_dist=0.0
    radius=0.5
    num_closest=5
    OutRedPkl=None
    # 参数读取
    for i in range(len(sys.argv)):
        # For KmerFrequencyStatistics()
        if   '-fa'     == sys.argv[i]:
            Fasta_Path  = sys.argv[i+1]  # 必填
        elif '-kmin'   == sys.argv[i]:
            Min_ngram   = int(sys.argv[i+1])  # 建议默认
        elif '-kmax'   == sys.argv[i]:
            Max_ngram   = int(sys.argv[i+1])  # 建议默认
        # For DimensionalReductionVectors()
        elif '-size'   == sys.argv[i]:
            SampleSize  = int(sys.argv[i+1])  # 建议填
        elif '-pre'    == sys.argv[i]:
            preprocess  = sys.argv[i+1]  # 建议默认
        elif '-mod'    == sys.argv[i]:
            method      = sys.argv[i+1]  # 建议默认
        elif '-seed'   == sys.argv[i]:
            random_state= int(sys.argv[i+1])  # 建议填 
        elif '-pkl'    == sys.argv[i]:
            OutRedPkl   = sys.argv[i+1]
        # For GetRepSeq()
        elif '-out'    == sys.argv[i]:
            Otpt_RepSeq = sys.argv[i+1]  # 必填
        elif '-otnum'  == sys.argv[i]:
            num_closest = int(sys.argv[i+1])  # 建议默认
        elif '-radius'  == sys.argv[i]:
            radius       = float(sys.argv[i+1]) 
    
    
    
    S=time.time()
    # 程序运行
    print('='*5+'Step-1: Vectorization of DNA sequence features'+'='*5)
    NameLst,SeqArr=KmerFrequencyStatistics(Fasta_Path,Min_ngram,Max_ngram)
    print('Done!')
    
    print('='*5+'Step-2: Data dimensionality reduction'+'='*5)
    ReducedArr=DimensionalReductionVectors(SeqArr,preprocess,method,n_components,random_state,min_dist)
    if OutRedPkl != None:
        print(f'Output Reduced_Array to {OutRedPkl}')
        with open(OutRedPkl,'wb') as f:
            pickle.dump((NameLst,ReducedArr),f)    
    print('Done!')
    
    print('='*5+'Step-3: Clustering and Representative Sequence Extraction'+'='*5)
    GetRepSeq(NameLst,ReducedArr,SampleSize,Fasta_Path,Otpt_RepSeq,radius,num_closest)
    print('Done!')
    
    E=time.time()
    timestamp = E-S
    total_seconds = int(timestamp)
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    print(f'TIME COST: {str(hours)}h {str(minutes)}m {str(seconds)}s')