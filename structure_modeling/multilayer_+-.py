#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2020/10/29
# @Function: build structure with specified stacking offset direction

from ase.io import read,write
import numpy as np
from ase import Atom
import pandas as pd
from ase.visualize import view
import os

def read_1st_layer(file_name, diff_1_2):
    a_0 = read(file_name)
    a = a_0[0:atom_num]
    atom_symbol = []
    atom_positions = []
    # read atoms in the 1st layer
    for n in range(len(a_0)):
        if a_0.get_positions()[n, 2] > diff_1_2:
            atom_symbol.append(a_0.get_chemical_symbols()[n])
            atom_positions.append(a_0.get_positions()[n])
    a.set_chemical_symbols(atom_symbol)
    a.set_positions(atom_positions)
    return a
def get_angle(a,b,c):
    La = np.sqrt(a.dot(a))
    Lb = np.sqrt(b.dot(b))
    Lc = np.sqrt(c.dot(c))
    cos_alpha = b.dot(c)/(Lb*Lc)
    cos_beta = a.dot(c)/(La*Lc)
    alpha = np.degrees(np.arccos(cos_alpha))
    beta = np.degrees(np.arccos(cos_beta))
    return alpha, beta
def control_c(vector, vector_name, c):
    if vector_name == 'a_vector':
        n = 0
    elif vector_name == 'b_vector':
        n = 1
    if c[n] < 0 and vector[n] < 0:
        c -= vector
    elif c[n] < 0 and vector[n] > 0:
        c += vector
    else:
        print("c_n is along vector a and b")
    return c
def def_cell(periodicity='bulk'):
    a1 = a.get_positions()[0]
    a2 = a.get_positions()[atom_num ]
    aN = a.get_positions()[atom_num * (N - 1)]
    d = a.get_positions()[atom_num][2] - a.get_positions()[0][2]
    a_vector = lattice[0]
    b_vector = lattice[1]
    c_n = aN - a1 + a2 - a1
    #c_n[2] += d
    return c_n
''' # check the angle of alpha and beta
    alpha = get_angle(a_vector, b_vector, c_n)[0]
    beta = get_angle(a_vector, b_vector, c_n)[1]
    print('before control_direction: ', alpha, beta, c_n)  # 41.44178862173217 90.53903320550464 [-0.15 15.79  9.35]

    c_n = control_c(a_vector, 'a_vector', c_n)
    c_n = control_c(b_vector, 'b_vector', c_n)
    print('after control_direction: ', alpha, beta, c_n)

    for i in range(1): # i=0,a_vector; i=1, b_vector. a is along x axis, b is along y axis
        ref = a.get_cell()[i]
        if ref[i] > 0 and c_n[i] > 0:
            if c_n[i] - ref[i] > ref[i]:
                c_n[i] = c_n[i] - ref[i] * int((c_n[i] - ref[i])/ref[i])
            else:
                print ('ref[i] > 0 and comp[i] > 0')
        elif ref[i] < 0 and c_n[i] < 0:
            if c_n[i] - ref[i] < ref[i]:
                c_n[i] = c_n[i] - ref[i] * int((c_n[i] - ref[i])/ref[i])
            else:
                print ('ref[i] < 0 and comp[i] < 0')
        elif ref[i] > 0 and c_n[i] < 0:
            if ref[i] - c_n[i] > ref[i]:
                c_n[i] = c_n[i] + ref[i] * int((ref[i] - c_n[i])/ref[i])
            else:
                print('ref[i] > 0 and comp[i] < 0')
        elif ref[i] < 0 and c_n[i] > 0:
            if c_n[i] - ref[i] > -1 * ref[i]:
                c_n[i] = c_n[i] + ref[i] * int((c_n[i] - ref[i])/ref[i])
            else:
                print('ref[i] < 0 and comp[i] > 0')
        else:
            print("can't comp", i, ref)

    alpha = get_angle(a_vector, b_vector, c_n)[0]
    beta = get_angle(a_vector, b_vector, c_n)[1]
    if periodicity == 'slab':
        c_n[2] += 15
    print('after control_angle: ', alpha,beta, c_n)
    return c_n'''


############################ Input parameters #################################

path = 'E:/1.1 stacking_COF/some_files/'
sta_data = pd.read_csv(path+'sta_list.csv')
start = 10 # 4l: [10:13]; 6l: [0:10]; 100l_inc & sta: [13:15]
end = 13
data = pd.read_csv(path+'Energies_COF-1.csv') # read energies
data_vector = pd.read_csv(path + 's_vector_COF-1.csv')
sta_list = sta_data['sta_seq']
file_name = 'AA_node_COF1.gen'
file_type = 'gen'
out_path = 'D:/project/COF-1/'

diff_1_2 = 2 # (for COF-1, AA: 2, AB: 0.1; for COF-5, 3; for graphyne: 0)
atom_num = 42

shift_types = ['eclipsed', 'node', 'linker', 'linker_m']
s_plus  = data_vector['AA_n'][0:3] # COF-1:[-0.0340,1.4810,3.5530]
s_minus = (s_plus[0], -s_plus[1], s_plus[2]) # COF-1:[0.0340,-1.4810,3.5530]
s_lattice = [0,0,0]

lattice = read(file_name).get_cell()
lattice_b = lattice[0,1]+lattice[1,1]
N = len(sta_list[start]) + 1
os.mkdir(out_path + str(N) + 'la_+-')

for j in range(start, end):
    #a = read(file_name)
    a = read_1st_layer(file_name, diff_1_2)
    a.bak = a[:]
    m = [0, 0, 0]
    b = []
    for i in range(atom_num):
        m[0] = a.get_positions()[i, 0]
        m[1] = a.get_positions()[i, 1]
        m[2] = a.get_positions()[i, 2] + 1
        m_ = m[:]
        b.append(m_)
    a.set_positions(b)
    #a = a.bak[:]
    k = 0

    sta_seq = sta_list[j]
    #sta_seq = sta_list[:]
    for i in sta_seq:
        if i == '+':
            s = s_plus
        elif i == '-':
            s = s_minus
        else:
            print(i,'error with sta_seq')
        #print(i,s)
        m = np.array([0, 0, 0])
        #s = s *2
        for i in range(atom_num * k, atom_num * (k + 1)):
            b = a.get_positions()
            a.append(Atom(a.get_chemical_symbols()[i]))
            m = a.get_positions()[i] + s_lattice + s
            a.set_positions(np.vstack((b, m)))
        k += 1
        #print(j, len(a))
    # set the cell for the N-layer model
    c_n = def_cell()
    a.set_cell([a.get_cell()[0], a.get_cell()[1], c_n])
    write(out_path + str(N) + 'la_+-' + '/'+ str(N) + '_' + str(sta_seq) + '.' + file_type, a, format=file_type)  # output the structure
#str(out_path)+ str(N) + '_' + str(sta_seq)+'/'
#p = read(out_filename, a, format=file_type)
#p=read(str(out_path)+ str(N) + '_' +  str(sta_seq)+'/'+'dftb_in.gen')
#view(p)

