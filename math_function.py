'''
Author: qyp422
Date: 2023-03-13 16:03:49
Email: qyp422@qq.com
LastEditors: Please set LastEditors
LastEditTime: 2023-03-13 20:02:16
Description: 

Copyright (c) 2023 by qyp422, All Rights Reserved. 
'''
import numpy as np
from numba import njit
#并查集
class Quick_Find():
    def __init__(self, n):
        self.count = n
        self._parent = [i for i in range(n)]
        self._weight = [1 for i in range(n)]

    def union(self, p, q):
        rootP = self.find(p)
        rootQ = self.find(q)
        if rootP == rootQ:
            return
        # 轻根到重根，为了平衡
        if self._weight[rootP] > self._weight[rootQ]:
            self._parent[rootQ] = rootP
            self._weight[rootP] += self._weight[rootQ]
        else:
            self._parent[rootP] = rootQ
            self._weight[rootQ] += self._weight[rootP]
        self.count -= 1

    def is_connected(self, p, q):
        return self.find(p) == self.find(q)

    def find(self, x):
        while self._parent[x] != x:
            # 路径压缩
            self._parent[x] = self._parent[self._parent[x]]
            x = self._parent[x]
        return x

    def get_count(self):
        return self.count
    
    
    '''
    description: 
    param {*} self
    return {*} result[0] 最大的cluster对应的id
    '''
    def hash_dict(self):
        hash1 = {}
        for i in range(len(self._parent)):
            parent = self.find(i)
            if parent in hash1:
                hash1[parent].append(i)
            else:
                hash1[parent] = [i]
        sort_id = []
        keys = []
        for j in hash1:
            sort_id.append(j)
            keys.append(len(hash1[j]))
        keys_sort =np.argsort(np.array(keys))
        result = [hash1[sort_id[x]] for x in keys_sort]
        return result.reverse()



'''
description: 所有的粒子到达某一粒子的距离一定要小于0.5box长度
            若这些粒子的最左边到最右边超过0.5box长度计算结果将会不准确
param {*} x x坐标np.array
param {*} y y坐标np.array
param {*} z z坐标np.array
param {*} mass mass坐标np.array
param {*} box_array
return {*}
'''
def cal_cm_rg(x,y,z,mass,box):
    n = len(mass)
    cm0 = np.array([x[0],y[0],z[0]]) #幸运粒子
    l_box = np.array([box[1] - box[0],box[3] - box[2],box[5] - box[4]],dtype=np.double)
    cm = np.zeros(3,dtype=np.double)
    pos = np.zeros((n,3),dtype=np.double)
    total_mass = np.sum(mass)
    for i in range(n):
        pos[i] = np.array([x[i],y[i],z[i]])
        dr = pos[i]-cm0
        pos[i] -= l_box * np.rint (dr / l_box) #pos 距离cm0为盒子一半以内
        cm += pos[i]*mass[i]
    cm /= total_mass

    rg_sq = 0.0
    for i in range(n):
        rg_sq += np.dot(pos[i] - cm,pos[i] - cm)*mass[i]
    rg_sq /= total_mass

    return cm,np.sqrt(rg_sq)

