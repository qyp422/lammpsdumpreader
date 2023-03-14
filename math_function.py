'''
Author: qyp422
Date: 2023-03-13 16:03:49
Email: qyp422@qq.com
LastEditors: Please set LastEditors
LastEditTime: 2023-03-14 14:24:25
Description: 

Copyright (c) 2023 by qyp422, All Rights Reserved. 
'''
import numpy as np
from numba import njit,vectorize

import logging
numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)
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
        result.reverse()
        return result



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


@njit
def cal_cm_rg(x,y,z,l_box,mass):
    n = len(x)
    cm0 = np.array([x[0],y[0],z[0]]) #幸运粒子
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

@vectorize
def dis_sq(a,b,l_box):
    dx = min(abs(a[0] - b[0]),l_box[0]-abs(a[0] - b[0]))
    dy = min(abs(a[1] - b[1]),l_box[1]-abs(a[1] - b[1]))
    dz = min(abs(a[2] - b[2]),l_box[2]-abs(a[2] - b[2]))
    return dx**2+dy**2+dz**2

'''
description: 
param {*} s numpy定义的“结构体”
param {*} l_box 盒子尺寸
param {*} num_chain 链条数目
param {*} num_beads 每条链的珠子数目
return {*}
'''
@njit
def get_mol_cm_rg(s,l_box,num_chain,num_beads):
    cm = np.zeros((num_chain,3),dtype=np.double)
    rg = np.zeros((num_chain,),dtype=np.double)
    for j in range(num_chain):
        cm[j],rg[j] = cal_cm_rg(s['x'][j*num_beads:j*num_beads+num_beads],s['y'][j*num_beads:j*num_beads+num_beads],s['z'][j*num_beads:j*num_beads+num_beads],l_box,s['mass'][j*num_beads:j*num_beads+num_beads])
    return cm,rg

def find_cluster(pair,num_chain):
    q = Quick_Find(num_chain)
    for (i,j) in pair:
        q.union(i,j)
    cluster_list = q.hash_dict()
    del q
    return cluster_list


def cluster_to_ids(mol_list,num_beads,start=0):
    res = []
    for i in mol_list:
        res += list(range(start+i*num_beads,start+i*num_beads+num_beads))
    return res


@njit
def cluster_ysz(cm,num_chain,cut_off,l_box):
    pair = {}
    for i in range(num_chain-1):
            for j in range(i+1,num_chain):
                dx = min(abs(cm[i][0] - cm[j][0]),l_box[0]-abs(cm[i][0] - cm[j][0]))
                dy = min(abs(cm[i][1] - cm[j][1]),l_box[1]-abs(cm[i][1] - cm[j][1]))
                dz = min(abs(cm[i][2] - cm[j][2]),l_box[2]-abs(cm[i][2] - cm[j][2]))
                if dx*dx + dy*dy + dz*dz <= cut_off*cut_off:
                    pair[(i,j)] = 1
    return pair

@njit
def cluster_sy(s,num_chain,num_beads,cut_off,l_box,type1,type2):
    pair = {}
    intrapair = np.array([0 for _ in range(num_chain)])
    outerpair = 0
    for i in type1:
        for j in type2:
            dx = min(abs(s['x'][i] - s['x'][j]),l_box[0]-abs(s['x'][i] - s['x'][j]))
            dy = min(abs(s['y'][i] - s['y'][j]),l_box[1]-abs(s['y'][i] - s['y'][j]))
            dz = min(abs(s['z'][i] - s['z'][j]),l_box[2]-abs(s['z'][i] - s['z'][j]))
            if dx*dx + dy*dy + dz*dz <= cut_off*cut_off:
                mol1 = i//num_beads
                mol2 = j//num_beads
                if mol1 == mol2:
                    intrapair[mol1] += 1
                else: 
                    if mol1 > mol2:
                        mol1,mol2 = mol2,mol1
                    if (mol1,mol2) in pair:
                        pair[(mol1,mol2)] += 1
                    else:
                        pair[(mol1,mol2)] = 1
                    outerpair += 1
                break

    return pair,intrapair,outerpair
