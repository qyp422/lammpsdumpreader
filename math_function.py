'''
Author: qyp422
Date: 2023-03-13 16:03:49
Email: qyp422@qq.com
LastEditors: Please set LastEditors
LastEditTime: 2023-03-13 16:46:57
Description: 

Copyright (c) 2023 by qyp422, All Rights Reserved. 
'''
import numpy as np

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

