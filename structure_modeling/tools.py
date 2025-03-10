#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2021/11/24
# @Function: small tools often use for structure modeling/transform/change

from ase.io import read,write
from ase.visualize import view
import os.path
import numpy as np

def conv2cif(file_name, path='structures/',ftype='cif'):
    fname, ext = os.path.splitext(file_name)
    a = read(path+file_name)
    out_file = fname + '.' + ftype
    write(path+out_file, a, format=ftype)
    return out_file
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
def def_cell():
    a1 = a.get_positions()[0]
    a2 = a.get_positions()[atom_num * (N - 1)]
    d = a.get_positions()[atom_num][2] - a.get_positions()[0][2]
    a_vector = lattice[0]
    b_vector = lattice[1]
    c_n = a2 - a1
    c_n[2] += d

    # check the angle of alpha and beta
    alpha = get_angle(a_vector, b_vector, c_n)[0]
    beta = get_angle(a_vector, b_vector, c_n)[1]
    print(alpha,beta, c_n)  # 41.44178862173217 90.53903320550464 [-0.15 15.79  9.35]

    c_n = control_c(a_vector, 'a_vector', c_n)
    c_n = control_c(b_vector, 'b_vector', c_n)
    print('after control_c', c_n)

    for i in range(1): # i=0,a_vector; i=1, b_vector. a is along x axis, b is along y axis
        ref = a.get_cell()[i]
        comp = c_n
        if ref[i] > 0 and comp[i] > 0:
            if comp[i] - ref[i] > ref[i]:
                c_n[i] = c_n[i] - ref[i]
            else:
                print ('ref[i] > 0 and comp[i] > 0')
        elif ref[i] < 0 and comp[i] < 0:
            if comp[i] - ref[i] < ref[i]:
                c_n[i] = c_n[i] - ref[i]
            else:
                print ('ref[i] < 0 and comp[i] < 0')
        elif ref[i] > 0 and comp[i] < 0:
            if comp[i] - ref[i] < -1 * ref[i]:
                c_n[i] = c_n[i] + ref[i]
            else:
                print('ref[i] > 0 and comp[i] < 0')
        elif ref[i] < 0 and comp[i] > 0:
            if comp[i] - ref[i] > -1 * ref[i]:
                c_n[i] = c_n[i] + ref[i]
            else:
                print('ref[i] < 0 and comp[i] > 0')
        else:
            print("can't comp", i, ref, comp)

    alpha = get_angle(a_vector, b_vector, c_n)[0]
    beta = get_angle(a_vector, b_vector, c_n)[1]
    print(alpha,beta, c_n)
    return c_n
def view_str(file_name, path='structures/'):
    p = read(path + file_name)
    view(p)


#path = "C:/Users/Qianqian_Zhang/Desktop/COF-1/pXRD_generation/ptept/AB_ABC_only_10_0.xyz"
#path = "D:/project/graphyne/2layer_bulk/AA_l_inclined/dftb+_UFF/"
path = "D:/project/COF-1/3layer_slab/AB_node/240_node/"

#out_path = 'D:/project/COF-5/2layer_slab/'
file_name = 'geo_end.gen'
#conv2cif(file_name,path=path)
#view_str('dftb_in.gen',path=out_path + 'AB_l' + '/')
view_str(file_name,path=path)


'''def set_cell(periodicity='bulk'):
    atom1 = a.get_positions()[0]
    atom2 = a.get_positions()[atom_num * (N - 1)]
    addi_space = a.get_positions()[atom_num]-a.get_positions()[0]
    a_vector = lattice[0]
    b_vector = lattice[1]
    c_n = atom2 - atom1 + addi_space
    #print('c_n', c_n)
    #print('a_vector',a_vector)
    #print('b_vector',b_vector)

# check
    if c_n[1] >= 0:
        while np.abs(c_n[1]) > np.abs(b_vector[1]):
            c_n -= b_vector
            #print('c_n[1]>0',c_n)
    elif c_n[1] < 0:
        while np.abs(c_n[1]) > np.abs(b_vector[1]):
            c_n += b_vector
            #print('c_n[1]<0',c_n)
    if c_n[0] >= 0:
        while np.abs(c_n[0]) > np.abs(a_vector[0]):
            c_n -= a_vector
            #print('c_n[0]>0',c_n)
    elif c_n[0] < 0:
        while np.abs(c_n[0]) > np.abs(a_vector[0]):
            c_n += a_vector
            #print('c_n[0]<0',c_n)

    # make c_n along +a and +b
    """if c_n[0] > 0 and a_vector[0] < 0:
        c_n[0] += a_vector[0]
    elif c_n[0] < 0 and a_vector[0] >0:
        c_n[0] += a_vector[0]

    if c_n[1] > 0 and a_vector[1] < 0:
        c_n[1] += a_vector[1]
    elif c_n[1] < 0 and a_vector[1] >0:
        c_n[1] += a_vector[1]"""

    # check the angle of alpha and beta
    alpha = np.arctan(c_n[2] / c_n[0])
    beta = np.arctan(c_n[2] / c_n[1])
    angle_ref = 85
    angle_ref = np.radians(angle_ref)
    #print(np.degrees(alpha), np.degrees(beta))
    if np.abs(alpha) < angle_ref:
        if np.abs(beta) < angle_ref:
            # make c_n along +a and +b
            if c_n[0] > 0 and a_vector[0] < 0:
                c_n[0] += a_vector[0]
            elif c_n[0] < 0 and a_vector[0] > 0:
                c_n[0] += a_vector[0]

            if c_n[1] > 0 and a_vector[1] < 0:
                c_n[1] += b_vector[1]
            elif c_n[1] < 0 and a_vector[1] > 0:
                c_n[1] += b_vector[1]
    if periodicity == 'slab':
        c_n[2] += 15
    return c_n
'''