# ForSklearn.py
import numpy as np

def calculate_centroid(points):
    centroid_x = np.median(points[:, 0])
    centroid_y = np.median(points[:, 1])
    return np.array([centroid_x, centroid_y])

def Get_DBSCAN_CenterPosition_00000(data,dbscan_rst):
    '''
    输入：DBSCAN的fit(array)结果
    返回：DBSCAN结果的各簇密度中心(密度最高的位置)
    '''
    labels = [x+1 for x in dbscan_rst.labels_]
    Clusters = sorted(list(set(labels))) # 有几个簇按从小到大排列
    max_density_points = {}
    for cluster in Clusters:
        if cluster != 0:  # 忽略噪声点（标记为0），原始结果应该是-1但我在上面+1了
            cluster_points = data[labels == cluster]   # 提取改簇的所有点的坐标数据
            # 计算簇内点的密度（可以根据需要选择其他距离度量方法）
            density = len(cluster_points)/(np.sum((cluster_points[:, np.newaxis]-cluster_points)**2,axis=2)<dbscan_rst.eps**2).sum()
            # 找到密度最大的坐标
            max_density_points[cluster]=cluster_points[np.argmax(density)]
    return max_density_points

def Get_DBSCAN_CenterPosition(data,dbscan_rst):
    '''
    输入：DBSCAN的fit(array)结果
    返回：DBSCAN结果的各簇密度中心(密度最高的位置)
    '''
    labels = [x+1 for x in dbscan_rst.labels_]
    Clusters = sorted(list(set(labels))) # 有几个簇按从小到大排列
    centers = {}
    for cluster in Clusters:
        if cluster != 0:  # 忽略噪声点（标记为0），原始结果应该是-1但我在上面+1了
            cluster_points = data[labels == cluster]   # 提取该簇的所有点的坐标数据
            total_num=len(cluster_points)    # 该簇所有点的数量
            centers[cluster]={'POS':calculate_centroid(cluster_points),'NUM':total_num}
    return centers