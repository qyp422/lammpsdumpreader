# lammpsdumpreader
lammps dumper文件读取与cluster分析
## dumpreader.py 
### class Lammps_dumpreader 用来读取lammpstrj文件 

关于质量的输入如果lammpstrj里面有mass一项则优先使用lammpstrj中的mass
其次会根据Lammps_dumpreader的m_array参数设置默认值1.0或者自己根据原子类型设置mass。

可选参数m_array若为float值则设置每个原子为固定的质量，若为dict或list需要确保m[typeid]为对应typeid的质量。默认值为1.0 ！
## analysis.py
用来执行分析脚本 请根据自己体系设置cut_off num_chain num_beads

使用方法为 python analysis.py lammpstrj_filename -bm -cm -r start end step -o output_filename -yzs/-sy -rg --c -mc -mcb -cv
### cluster方案 -yzs/-sy
对应不同的体系找cluster方案，也可以后期diy

### 可选参数 -bm 
将盒子质心移动至[0,0,0]，并没有移动粒子！若想粒子也平移请使用shift_pos[0,0,0]!
### 可选参数 -cm
将最大的cluster质心移动至盒子质心

### 可选参数 -r 
不设置即读取所有帧

若为1个参数即读取帧数为第0帧到end(包括end) eg -r end

若为2个参数即读取帧数为start帧到end(包括end) eg -r start end

若为3个参数即读取帧数为start帧到end，每隔step读取一次(包括end) eg -r start end step

### 可选参数 -rg 
输出每个分子的rg 至文件名lammpstrj_filename_rg.txt中，格式见输出文件
### 可选参数 -c
输出体系所有非单体cluster 至文件名cluster.txt中，格式见输出文件
输出-c中的cluster.txt统计结果至文件名drop.txt中，格式见输出文件
### 可选参数 -mc
输出体系最大的cluster 至文件名max_cluster.txt中，格式见输出文件
### 可选参数 -cv
统计每个cluster大小zise从1-9，>=10至文件名cluster_volum.txt中，格式见输出文件
### 可选参数 -mcb
输出体系最大的cluster的成键信息，至文件名max_cluster_bonds.txt中，格式见输出文件

注意必须在-sy体系下
### 可选参数 -o
-o output_filename

输出分析的轨迹，比如做过-bm -cm -r 将会输出用作分析后的轨迹文件,新的轨迹文件名为output_filename。

## math_function.py
数学计算公式
### class Quick_Find
快速并查集并可以按照cluster从大到小输出分组
### def cal_cm_rg
计算周期性边界条件下的质心坐标已经回转半径rg
### dis_sq
计算两个向量在周期性边界条件下的距离平方
### def get_mol_cm_rg(s,l_box,num_chain,num_beads)
寻找分子的cm与rg，返回np.array

## 输出文件 
### bond_number.txt
timestep 分子外成键数 分子内成键数 成键总数 未成键位点
### lammpstrj_filename_rg.txt
timestep 第一个分子rg 第二个分子rg …………
### cluster.txt
timestep = timestep

cluster1_rg cluster1大小 （cluster1组成的分子id-1）

cluster2_rg cluster2大小 （cluster2组成的分子id-1）

cluster3_rg cluster3大小 （cluster3组成的分子id-1）

按照cluster大小从大到小排列
### max_cluster.txt
timestep cluster大小 （组成的分子id-1）
### max_cluster_bonds.txt
timestep 分子外成键数 分子内成键数 成键总数 未成键位点

note只统计最大的cluster
### drop.txt
timestep 平均cluster_rg cluster个数 cluster大小之和  
### cluster_volum.txt
timestep 单聚体个数 2聚体个数 ……………… 9聚体个数 超过9聚体个数 所有超过9聚体链条数总和

