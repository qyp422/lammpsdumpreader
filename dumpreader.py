'''
Author: qyp422
Date: 2022-10-17 15:29:44
Email: qyp422@qq.com
LastEditors: Please set LastEditors
LastEditTime: 2023-03-13 20:10:06

Description: 

Copyright (c) 2022 by qyp422, All Rights Reserved. 
'''


import numpy as np

print('\n' + __file__ + ' is called\n')

class Lammps_dumpreader():
    def __init__(self, configuration,m_array = 1.):
        self.frame_num = 0
        self._timestep = -1 
        self._conf = False
        self.total_atoms = False
        self.box_arr = False
        self.l_box = False
        try:
            self._conf = open(configuration, "r")
        except:
            exit(f'cannot find file {configuration}')

        self._system = False
        self._component = ['mol','type','q','x','y','z','mass']
        self._dumpid = {}
        self._m_array = m_array
        self._type_index = {}

        if self._update_system_parameter():
            print('Lammps_dumpreader init done!!!!!\n')
            self._conf.close() 
            self._conf = open(configuration, "r")
        else:
            self._conf.close() 
            exit('dumpfile formate is wront!')

    def __del__(self):
        if self._conf:self._conf.close()
    
    #读取一帧，可以跳过数据
    def _read_single_frame(self, skip=False):
        line = self._conf.readline()
        if  line == '':
            print('\nend of data reading!\n')
            return False
        if skip:
            # need to skip elf.total_atoms+8 lines
            self._timestep = int(self._conf.readline())
            for _ in range(self.total_atoms+7):
                self._conf.readline()
            print(f'data skip at system_time {self._timestep},total frame is {self.frame_num}!')
            self.frame_num += 1
            return True
        #   ITEM: TIMESTEP
        self._timestep = int(self._conf.readline())
        #   ITEM: NUMBER OF ATOMS
        self._conf.readline()
        # total_atoms
        self._conf.readline()
        #   ITEM: box BOUNDS pp pp pp
        self._conf.readline()
        # box xyz
        self._conf.readline()
        self._conf.readline()
        self._conf.readline()

        
        #   ITEM: ATOMS id mol type x y z c_quat[1] c_quat[2] c_quat[3] c_quat[4]
        self._conf.readline()
        for i in range(self.total_atoms):
            ls = self._conf.readline().split()
            self._system['x'][i] = float(ls[self._dumpid['x']])
            self._system['y'][i] = float(ls[self._dumpid['y']])
            self._system['z'][i] = float(ls[self._dumpid['z']])
        
        print(f'data updata at system_time {self._timestep},total frame is {self.frame_num}!')
        self.frame_num += 1
        # print(self._system)
        return True

    #初始化系统
    def _update_system_parameter(self):
        line = self._conf.readline()
        if  line == '':
            print('\nno data!\n')
            return False
        #   ITEM: TIMESTEP
        self._conf.readline()
        #   ITEM: NUMBER OF ATOMS
        self._conf.readline()
        self.total_atoms = int(self._conf.readline())
        #   ITEM: box BOUNDS pp pp pp
        line = self._conf.readline()
        self.box_arr  =  [float(x) for x in self._conf.readline().split()] #x
        self.box_arr  += [float(x) for x in self._conf.readline().split()] #y   
        self.box_arr  += [float(x) for x in self._conf.readline().split()] #z
        self.l_box = np.array([self.box_arr[1] - self.box_arr[0],self.box_arr[3] - self.box_arr[2],self.box_arr[5] - self.box_arr[4]],dtype=np.double)
        #   ITEM: ATOMS id mol type x y z c_quat[1] c_quat[2] c_quat[3] c_quat[4]
        line = [x for x in self._conf.readline().split()]
        for k in range(2,len(line)):    #0 1 is ITEM: ATOMS
            self._dumpid[line[k]] = k-2
        if not all(['x' in self._dumpid,'y' in self._dumpid,'z' in self._dumpid]):
            print('cannot find x y z in line ITEM: ATOMS!')
            return False
        particle_dtype = np.dtype({'names':['mol','type','q','x','y','z','mass'], 
                                    'formats':[np.int64, 
                                               np.int64, 
                                               np.double, 
                                               np.double, 
                                               np.double,
                                               np.double,
                                               np.double]})
        self._system = np.zeros((self.total_atoms), dtype=particle_dtype)
        for i in range(self.total_atoms):
            ls = self._conf.readline().split()
            if 'mol' in self._dumpid:
                self._system['mol'][i] = int(ls[self._dumpid['mol']])
            if 'type' in self._dumpid:
                self._system['type'][i] = int(ls[self._dumpid['type']])
                # make a list of every type for quick union
                if self._system['type'][i] in self._type_index:
                    self._type_index[self._system['type'][i]].append(i)
                else:
                    self._type_index[self._system['type'][i]] = [i]
            if 'q' in self._dumpid:
                self._system['q'][i] = float(ls[self._dumpid['q']])
            if 'mass' in self._dumpid:
                self._system['mass'][i] = float(ls[self._dumpid['mass']])
            elif not isinstance(self._m_array,float):
                self._system['mass'][i] = float(self._m_array[int(ls[self._dumpid['type']])])
            else:
                self._system['mass'][i] = self._m_array
        print(f'total_atoms is {self.total_atoms}')
        return True

    #将盒子质心移动至原点
    def box_mid(self):
        print(f'old box array is {self.box_arr}')
        self.box_arr = [-(self.box_arr[1] - self.box_arr[0])/2,(self.box_arr[1] - self.box_arr[0])/2,
                        -(self.box_arr[3] - self.box_arr[2])/2,(self.box_arr[3] - self.box_arr[2])/2,
                        -(self.box_arr[5] - self.box_arr[4])/2,(self.box_arr[5] - self.box_arr[4])/2]
        print(f'new box array is {self.box_arr}')
        print(' ')

    #返回特定帧数的数据
    def get_frame_system(self,frame_num):
        if frame_num == 0:
            self._read_single_frame(skip = False)
            return self._system
        while self._read_single_frame(skip = True):
            if self.frame_num == frame_num:
                self._read_single_frame(skip = False)
                return self._system
        print(f'cannot found {frame_num} th frame! ')
        return False

    #粒子平移shift距离向左为-向右为+ 并根据周期性边界条件放置盒子中
    def shift_pos(self,shift):
        self._system['x']  = self._system['x'] + shift[0] - (self.box_arr[1]+self.box_arr[0])/2
        self._system['x']  = self._system['x']  - np.rint(self._system['x'] /(self.box_arr[1]-self.box_arr[0]))*(self.box_arr[1]-self.box_arr[0]) + (self.box_arr[1]+self.box_arr[0])/2

        self._system['y']  = self._system['y'] + shift[1] - (self.box_arr[3]+self.box_arr[2])/2
        self._system['y']  = self._system['y']  - np.rint(self._system['y'] /(self.box_arr[3]-self.box_arr[2]))*(self.box_arr[3]-self.box_arr[2]) + (self.box_arr[3]+self.box_arr[2])/2

        self._system['z']  = self._system['z'] + shift[2] - (self.box_arr[5]+self.box_arr[4])/2
        self._system['z']  = self._system['z']  - np.rint(self._system['z'] /(self.box_arr[5]-self.box_arr[4]))*(self.box_arr[5]-self.box_arr[4]) + (self.box_arr[5]+self.box_arr[4])/2

    #输出lammpstrj文件
    def out_put_lammpstrj(self,opt = ['mol','type','x','y','z']):
        out = ''
        out += 'ITEM: TIMESTEP\n'
        out += f'{self._timestep}\n'
        out += 'ITEM: NUMBER OF ATOMS\n'
        out += f'{self.total_atoms}\n'
        out += 'ITEM: BOX BOUNDS pp pp pp\n'
        out += f'{self.box_arr[0]} {self.box_arr[1]}\n'
        out += f'{self.box_arr[2]} {self.box_arr[3]}\n'
        out += f'{self.box_arr[4]} {self.box_arr[5]}\n'
        out += 'ITEM: ATOMS id '
        out += ' '.join(opt)
        out += '\n'
        for i in range(self.total_atoms):
            out += f'{i+1} ' + ' '.join([str(self._system[x][i]) for x in opt]) + '\n'
        return out