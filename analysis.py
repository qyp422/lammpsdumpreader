'''
Author: qyp422
Date: 2023-03-13 12:41:00
Email: qyp422@qq.com
LastEditors: Please set LastEditors
LastEditTime: 2023-03-13 17:23:28
Description:  

Copyright (c) 2023 by qyp422, All Rights Reserved. 
'''
import os,sys,datetime
import dumpreader as dr
import argparse
import math_function as mf

parser = argparse.ArgumentParser()
parser.add_argument("lammpstrj", help="input full name of .lammpstrj file!")
parser.add_argument("-bm",'--box_mid', help="let box center [0,0,0]",action="store_true")
parser.add_argument("-cm",'--cluster_mid', help="let max cluster center of mass [0,0,0]",action="store_true")
parser.add_argument("-o",'--output', help="output new lammpstrj",type=str)
parser.add_argument("-r",'--read', help="read freame from start to end every step",action="extend", nargs="+", type=int)


args = parser.parse_args()
#配置-r
frame_strat = 0
frame_end = -1
frame_step = 1
if not args.read:
    print('note read all frame!')
elif len(args.read) == 3:
    frame_strat = args.read[0]
    frame_end = args.read[1]
    frame_step = args.read[2]
elif  len(args.read) == 2:
    frame_strat = args.read[0]
    frame_end = args.read[1]
elif  len(args.read) == 1:
    frame_end = args.read[0]
else:
    exit('wrong -r commend!')




def main():

    start=datetime.datetime.now()

    mass1=[-1,131.199997,71.0800018,87.0800018,114.099998,115.099998,163.199997,101.099998,128.100006,57.0499992,97.1200027]
    # open file
    r = dr.Lammps_dumpreader(args.lammpstrj,m_array=1.) # m_array=mass1
    if args.box_mid:
        r.box_mid()

    if frame_strat == 0:
        condition = True #frame 0 read
    else:
        condition = False

    while r._read_single_frame(skip = not condition):
    
        if condition:
            
            if args.cluster_mid:
                r.shift_pos([0,0,0])
            if args.output:
                output.write(r.out_put_lammpstrj())
        # reader condition
        
        condition = True  if (not args.read) or ((frame_strat<= r.frame_num <=frame_end) and ((r.frame_num - frame_strat) % frame_step == 0)) else False

    del r
    end=datetime.datetime.now()
    print(f'start at {start}')
    print(f'end   at {end}')
    print(f'total wall time is  {end-start}')

    
if __name__ == "__main__":

    #open file
    # -o
    if args.output:
        output = open(args.output,'w+')
    

    
    ##############################
    main()
    ##############################

    #close file 
    if args.output:
        output.close