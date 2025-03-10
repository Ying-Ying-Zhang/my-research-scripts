#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2020/08/13
# @Function: build 2-layer structure with all possible stackings

from ase.io import read,write
import numpy as np
from ase import Atom
import pandas as pd
from ase.visualize import view
import os
from pystr import view_str

def add_2nd(shift_x,shift_y,shift_z):
    f_bak = []
    for i in range(atom_num):
        f_bak = f.get_positions()
        f.append(Atom(f.get_chemical_symbols()[i]))
        m[0] = f.get_positions()[i, 0] + shift_x
        m[1] = f.get_positions()[i, 1] + shift_y
        m[2] = f.get_positions()[i, 2] + shift_z
        f.set_positions(np.vstack((f_bak,m)))
    f.set_cell([f.get_cell()[0], f.get_cell()[1], (0,0,20)])
    return f

atom_num = 96
m = [0, 0, 0]
outf_type = 'gen'
in_path = 'D:/project/COF-5/1layer/dftb+_UFF/'
out_path = 'D:/project/COF-5/2layer_slab/'
file_name = 'geo_end.gen'
stacking_list = ['AA_e','AA_n','AA_l','AB_e','AB_n','AB_l']

for i in stacking_list:
    f = read(in_path + file_name)
    if i == 'AA_e':
        output = add_2nd(0,0,3.4)
    elif i == 'AA_n':
        output = add_2nd(0, 1.4, 3.2) #1/3*f.get_cell()[0,0]
    elif i == 'AA_l':
        output = add_2nd(1, 1, 3.2)
    elif i == 'AB_e':
        output = add_2nd(0,0 + 2 * f.get_cell()[1,1] / 3, 3.4)
    elif i == 'AB_n':
        output = add_2nd(0, 1.4 + 2 * f.get_cell()[1,1] / 3, 3.2)
    elif i == 'AB_l':
        output = add_2nd(1, 1 + 2 * f.get_cell()[1,1] / 3, 3.2)
    else:
        print('error with ' + i)
    os.mkdir(out_path + str(i))
    write(out_path + str(i) + '/' + 'dftb_in.gen', f, format = outf_type)

    #view_str('dftb_in.gen',out_path + str(i) + '/')