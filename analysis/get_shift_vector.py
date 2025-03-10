#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2020/12/02
# @Function: calculate the shift vector between layers 

from ase import Atom
from ase.io import read
import pandas as pd
import os

def s_vector(num_atom1,num_atom2):
    atom1 = a.get_positions()[num_atom1]
    atom2 = a.get_positions()[num_atom2]
    s = atom2 - atom1
    return s

path = 'D:/project/COF-1/3layer_slab/AB_linker/'
file_name = 'geo_end.gen'
folder_name = []
s1_x = []
s1_y = []
s1_z = []
s2_x = []
s2_y = []
s2_z = []
for files_folders in os.walk(path):  #List all folders and files in the current directory, [0]: current directory, [1]: folder, [2]: file
    #print(files_folders)
    for folders in files_folders[1]:
        #print('folders are:', folders)
        if folders != 'band' and 'geo_str' and 'in_str' and 'out_str':
            for file in os.listdir(path+folders):
                #print('file is:', file)
                if file == 'geo_end.gen':
                    folder_name.append(folders)
                    #file_path = os.path.join(path+folders, file)
                    a = read(path+folders+'/'+file)
                    s1_x.append(s_vector(1, 43)[0])
                    s1_y.append(s_vector(1, 43)[1])
                    s1_z.append(s_vector(1, 43)[2])
                    s2_x.append(s_vector(43, 85)[0])
                    s2_y.append(s_vector(43, 85)[1])
                    s2_z.append(s_vector(43, 85)[2])
                    #print(folder_name,s1_x)
#print(folder_name)
#print(len(folder_name),len(s1_x))
data = {'folder':folder_name,'s1_x':s1_x,'s1_y':s1_y,'s1_z':s1_z,'s2_x':s2_x,'s2_y':s2_y,'s2_z':s2_z}
f = pd.DataFrame(data, columns=['folder','s1_x','s1_y','s1_z','s2_x','s2_y','s2_z'])
f.to_csv(path+'shift_vector.csv')
