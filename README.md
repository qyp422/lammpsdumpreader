# lammpsdumpreader
lammps dumper文件读取与分析
## dumpreader.py 
class Lammps_dumpreader 用来读取lammpstrj文件 

可选参数m_array若为float值则设置每个原子为固定的质量，若为dict或list需要确保m[typeid]为对应typeid的质量。默认值为1.0 ！
## analysis.py
用来执行分析脚本

使用方法为 python analysis.py

可选参数 -bm 将盒子质心移动至[0,0,0]，并没有移动粒子！
