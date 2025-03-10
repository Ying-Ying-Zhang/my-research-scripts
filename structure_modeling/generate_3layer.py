#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2020/02/12
# @Function: build 3-layer structure with all possible stackings

from ase.io import read,write
import numpy as np
import math
from ase import Atom
from ase.visualize import view
from scipy.optimize import fsolve
import pandas as pd
import os

# define the mirror function
def mirror_point(i):
    mirror_k = (mirror[1]/mirror[0])
    x,y = i[0],i[1]
    return [
            v_all[2,1]+y-mirror_k*(v_all[2,0]+x),
            mirror_k*(y-v_all[2,1])+x-v_all[2,0]
            ]
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
def def_cell(periodicity='bulk'):
    a1 = a.get_positions()[0]
    a2 = a.get_positions()[atom_num ]
    aN = a.get_positions()[atom_num * (N - 1)]
    #d = a.get_positions()[atom_num][2] - a.get_positions()[0][2]
    c_n = aN - a1 + a2 - a1
    #c_n[2] += d

    l_a = a.get_cell()[0]
    l_b = a.get_cell()[1]
    print('before',c_n,l_a,l_b)
    while np.abs(c_n[1])>0.8*np.abs(l_b[1]):
        if (c_n[1]>0 and l_b[1]>0) or (c_n[1]<0 and l_b[1]<0):
            c_n-=l_b
            #print('1')
        elif (c_n[1] < 0 and l_b[1] > 0) or (c_n[1] > 0 and l_b[1] < 0):
            c_n += l_b
            #print('2')
    while np.abs(c_n[0])>0.8*np.abs(l_a[0]):
        if (c_n[0]>0 and l_a[0]>0) or (c_n[0]<0 and l_a[0]<0):
            c_n-=l_a
            #print('3')
        elif (c_n[0] < 0 and l_a[0] > 0) or (c_n[0] > 0 and l_a[0] < 0):
            c_n += l_a
            #print('4')
    if c_n[1]<0 and np.arctan(np.abs(c_n[2]/c_n[1]))<np.radians(80):
        c_n+=l_b
    if c_n[0]<0 and np.arctan(np.abs(c_n[2]/c_n[0]))<np.radians(80):
        c_n+=l_a
    print('after',c_n)
    return c_n

############################ Input parameters #################################
N = 3
n_rot = 6
r = 0
diff_1_2 = 0 # for input file fo AA: 2, AB: 0.1

path = 'E:/1.1 stacking_COF/some_files/'
path_gen = "D:/project/COF-5/2layer_slab/AB_linker/dftb+_UFF/"
data_vector = pd.read_csv(path + 's_vector_COF-1.csv') # read dat
stacking_type = 'AB' # stacking_type for teh 2layer structure: AA, AB
shift_type = 'linker' # shift vector of the provided file: node or linker
file_name = 'geo_end.gen'
generate_type = 'gen'
shift_types = ['eclipsed', 'node', 'linker', 'linker_m']
atom_num = 42
a = read_1st_layer(path_gen+file_name, diff_1_2)

lattice = a.get_cell()
lattice_b = lattice[0,1]+lattice[1,1]
s_AB_lattice = [0, 2 * lattice_b / 3, 0]
s_AA_lattice = [0,0,0]
s_AA = np.array([data_vector['AA_e'][0:3],data_vector['AA_n'][0:3],data_vector['AA_l'][0:3]])  # shift along node and linker directions
mirror_AA = data_vector['AA_l_mirror']
s_AB = np.array([data_vector['AB_e'][0:3],data_vector['AB_n'][0:3],data_vector['AB_l'][0:3]])   # eclipsed, node, linker
mirror_AB = data_vector['AB_l_mirror']

v_all = s_AA
mirror = mirror_AA
solution = fsolve(mirror_point, [0, 0])
u1 = np.array([solution[0], solution[1], s_AA[2, 2]])
s_AA = [s_AA[0], s_AA[1], s_AA[2], u1]

v_all = s_AB
mirror = mirror_AB
solution = fsolve(mirror_point, [0, 0])
u2 = np.array([solution[0], solution[1], s_AB[2, 2]])
s_AB = [s_AB[0], s_AB[1], s_AB[2], u2]

v_all = s_AA + s_AB # 0-3: AA_e,AA_n,AA_l,AA_l_m; # 4-7: AB_e,AB_n,AB_l,AB_l_m

if stacking_type == 'AA':
    #v_all = np.array([[0.0000,0.0000,3.5620],[-0.0340,1.4810,3.5530],[0.6510,1.3320,3.4790]]) # shift along node and linker directions
    s1_lattice = s_AA_lattice
    v_all_para = 0
elif stacking_type == 'AB':
    #v_all = np.array([[0, 0, 3.42], [-0.0130,1.1740,3.3490],[0.9590,1.2470,3.3570]])  # shift along node and linker directions
    s1_lattice = s_AB_lattice
    v_all_para = 4

#################################################################################
m = [0,0,0]
b=[]
for i in range(atom_num):
    m[0] = a.get_positions()[i, 0]
    m[1] = a.get_positions()[i, 1]
    m[2] = a.get_positions()[i, 2] + 0.5
    m_=m[:]
    b.append(m_)
a.set_positions(b)

if shift_type == 'eclipsed':
    shift_type_num = 3
    v = v_all[0+v_all_para]
    n_rot = 1
elif shift_type == 'node':
    shift_type_num = 3  # there are node, linker for the third layer of AA_node
    v = v_all[1+v_all_para]  # the shift vector for the 2nd layer
elif shift_type == 'linker':
    shift_type_num = 4  # there are node, linker, linker_m for the third layer of AA_linker
    v = v_all[2+v_all_para]  # the shift vector for the 2nd layer
else:
    print('Please input the correct shift_type')

for l in range(atom_num): # generate second layer
    b = a.get_positions()
    a.append(Atom(a.get_chemical_symbols()[l]))
    m = a.get_positions()[l] + s1_lattice + v
    a.set_positions(np.vstack((b, m)))
a_bak = a[:]

# generate AA/AB
os.mkdir(path + '3layer_slab/' + stacking_type + '_' + shift_type)
out_path = path + '3layer_slab/' + stacking_type + '_' + shift_type
## AA
s_vector_rotate_x = []
s_vector_rotate_y = []
s_vector_rotate_z = []
s_rotate_degree = []
for j in range(shift_type_num): # j from 0 if include eclipsed, otherwise j is from 1
    v = v_all[j]
    #print(v)
    name = shift_types[j]
    #print(stacking_type, n_rot,name)
    for i in range(n_rot):
        a = a_bak[:]
        r = 2 * i * np.pi / 6
        rotate = ([math.cos(r),math.sin(r)],[-math.sin(r),math.cos(r)])
        v1 = np.dot(rotate,([v[0],v[1]])) # rotate clockwise
        m = [0,0,0]
        s_vector_rotate_x.append(np.around(v1[0],3))
        s_vector_rotate_y.append(np.around(v1[1],3))
        s_vector_rotate_z.append(np.around(v[2],3))
        s_rotate_degree.append('AA_'+name+'_'+str(np.around(np.degrees(r),3)))
        for l in range(atom_num,atom_num*2):
            b=a.get_positions()
            a.append(Atom(a.get_chemical_symbols()[l]))
            m[0]=a.get_positions()[l,0]+s_AA_lattice[0]+v1[0]
            m[1]=a.get_positions()[l,1]+s_AA_lattice[1]+v1[1]
            m[2]=a.get_positions()[l,2]+s_AA_lattice[2]+v[2]
            a.set_positions(np.vstack((b,m)))
        a1 = a.get_positions()[0]
        a2 = a.get_positions()[atom_num]
        a3 = a.get_positions()[atom_num * 2]
        c = def_cell()
        #c = a3 - a1
        c[2] = c[2] + 15 #(a2[2] - a1[2])
        a.set_cell([a.get_cell()[0],a.get_cell()[1],c])
        write(out_path + '/' + str(60*i)+'_AA_'+str(name)+'_'+str(stacking_type)
              +'_'+str(shift_type)+'.'+generate_type,a,format=generate_type)  # output the 3layer

## AB
for j in range(shift_type_num): # j from 0 if include eclipsed, otherwise j is from 1
    v = v_all[j+4]
    print('v in AB', v)
    name = shift_types[j]
    #print(stacking_type, n_rot,name)
    for i in range(n_rot):
        a = a_bak[:]
        r = 2 * i * np.pi / 6
        rotate = ([math.cos(r),math.sin(r)],[-math.sin(r),math.cos(r)])
        v1 = np.dot(rotate,([v[0],v[1]])) # rotate clockwise
        m = [0,0,0]
        s_vector_rotate_x.append(np.around(v1[0],3))
        s_vector_rotate_y.append(np.around(v1[1]-s_AB_lattice[1],3))
        s_vector_rotate_z.append(np.around(v[2],3))
        s_rotate_degree.append('AB_'+name+'_'+str(np.around(np.degrees(r),3)))
        for l in range(atom_num,atom_num*2):
            b=a.get_positions()
            a.append(Atom(a.get_chemical_symbols()[l]))
            m[0]=a.get_positions()[l,0]-s_AB_lattice[0]+v1[0]
            m[1]=a.get_positions()[l,1]-s_AB_lattice[1]+v1[1]
            m[2]=a.get_positions()[l,2]-s_AB_lattice[2]+v[2]
            a.set_positions(np.vstack((b,m)))
        a1 = a.get_positions()[0]
        a2 = a.get_positions()[atom_num]
        a3 = a.get_positions()[atom_num * 2]
        c = def_cell()
        #c = a3 - a1
        c[2] = c[2] + 15 #(a2[2] - a1[2])
        a.set_cell([a.get_cell()[0],a.get_cell()[1],c])
        write(out_path + '/' + str(60*i)+'_AB_'+str(name)+'_'+str(stacking_type)
              +'_'+str(shift_type)+'.'+generate_type,a,format=generate_type)  # output the 3layer

# generate ABC
if stacking_type == 'AB':
    stacking_type = 'ABC'
    for j in range(shift_type_num): # j from 0 if include eclipsed, otherwise j is from 1
        v = v_all[j+4]
        name = shift_types[j]
        #print(stacking_type, shift_type_num, n_rot,name)
        for i in range(n_rot):
            a = a_bak[:]
            r = 2 * i * np.pi / 6
            rotate = ([math.cos(r),math.sin(r)],[-math.sin(r),math.cos(r)])
            v1 = np.dot(rotate,([v[0],v[1]])) # rotate clockwise
            m = [0,0,0]
            s_vector_rotate_x.append(np.around(v1[0], 3))
            s_vector_rotate_y.append(np.around(v1[1]-s_AB_lattice[1], 3))
            s_vector_rotate_z.append(np.around(v[2], 3))
            s_rotate_degree.append('ABC_' + name+'_'+ str(np.around(np.degrees(r), 3)))
            for l in range(atom_num,atom_num*2): # generate third layer ABC
                b=a.get_positions()
                a.append(Atom(a.get_chemical_symbols()[l]))
                m[0]=a.get_positions()[l,0]+s_AB_lattice[0]+v1[0]
                m[1]=a.get_positions()[l,1]+s_AB_lattice[1]+v1[1]
                m[2]=a.get_positions()[l,2]+s_AB_lattice[2]+v[2]
                a.set_positions(np.vstack((b,m)))
            a1 = a.get_positions()[0]
            a2 = a.get_positions()[atom_num]
            a3 = a.get_positions()[atom_num * 2]
            c = def_cell()
            #c = a3 - a1
            c[2] = c[2] + 15 # (a2[2] - a1[2])
            print(c[2])
            a.set_cell([a.get_cell()[0],a.get_cell()[1],c])
            write(out_path + '/' + str(60 * i) + '_ABC_' + str(name) + '_' + str(stacking_type)
                  + '_' + str(shift_type) + '.' + generate_type, a, format=generate_type)  # output the 3layer
#print(s_vector_rotate_z)
data = {'s_rotate_degree':s_rotate_degree,'s_vector_rotate_x':s_vector_rotate_x,'s_vector_rotate_y':s_vector_rotate_y,'s_vector_rotate_z':s_vector_rotate_z}
f = pd.DataFrame(data,columns=['s_rotate_degree','s_vector_rotate_x','s_vector_rotate_y','s_vector_rotate_z'])
f.to_csv(out_path+'_generate_shift_vector.csv')
#p = read(out_path + '/' + str(60 * i) + '_' + str(name) + '_' + str(stacking_type)
                  #+ '_' + str(shift_type) + '.' + generate_type)
#view(p)
