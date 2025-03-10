#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2024/12/28
# @Function: plot band structure and dos from aims output

import numpy as np
import pandas as pd
import re
import os
from math import modf
import fnmatch
import matplotlib.pyplot as plt


def set_xtick_dos(nrow):
    if 0 <= nrow < 10:
        xtick_dos = 5
    elif 10 <= nrow < 30:
        xtick_dos = 10
    elif 30 <= nrow < 60:
        xtick_dos = 20
    elif 60 <= nrow < 90:
        xtick_dos = 30
    elif 90 <= nrow < 150:
        xtick_dos = 50
    elif 150 <= nrow < 300:
        xtick_dos = 100
    else:
        xtick_dos = 200
    return xtick_dos


path = input('path is:') #'D:/project/COF-1/6l_bulk/6_+-+-+-/aims_pbe0/'
path_type = 'i'  # input('stacking is m, s i: ')   # monolayer; staggered; incliend
shift = '-'  # input('shift is + or -: ')
print(path)

add_num = 100
E_f = -4.87052550  # os.popen("grep 'Fermi level' aims.out | tail -1 |  awk '{print $6}' ")
Gap = 3.33760973  # os.popen("grep 'HOMO-LUMO gap' aims.out | tail -1 | awk '{print $5}' ")
e_f = float(E_f)
gap = float(Gap)
print('Fermi level is at: ', e_f)
print('band gap is: ', gap)

if shift == '+':
    E_shift = gap / 2
    print('s', 'E_shift is:', E_shift)
elif shift == '-':
    E_shift = -gap / 2
    print('i or m', 'E_shift is:', E_shift)

E_min = -round(gap/2 + 1.5, 0)
E_max = round(gap/2 + 1.5, 0)
print('E_min and E_max are: ', E_min, E_max)
if os.path.exists(path + '/band_aims' + str(round(gap,2)) + '.dat'):
    os.remove(path + '/band_aims' + str(round(gap,2)) + '.dat')
if os.path.exists(path + '/dos.dat'):
    os.remove(path + '/dos.dat')

band_files = fnmatch.filter(os.listdir(path),'band1*.out')
file_nums = []
for i in band_files:
    file_nums.append(re.findall(r'\d',i)[3])
file_num = max(file_nums)
print(band_files)
f1 = open(path + '/remain.dat','w+')
f2 = open(path + '/delete.dat','w+')
#-------------------------------get band data----------------------------------
for n in range(1,int(file_num)+1):
    file = 'band100' + str(n) + '.out'
    f = open(path + '/' + file)
    print(file)
    print(file + '\n', file = f1)
    print(file + '\n', file = f2)
    lines = f.readlines()
    raws_num = len(lines)  # 行数
    cols_num = len(re.findall(r'.[0-9]*[.][0-9]*', lines[0]))  # 列数

    band_data = np.zeros([raws_num, cols_num])
    for i in range(raws_num):
        line = re.findall(r'.[0-9]*[.0-9][0-9]*', lines[i])  # .[0-9]*[.0-9][0-9]*
        for j in range(cols_num):
            band_data[i, j] = float(line[j]) + float(E_shift)
    rows_num_data = np.shape(band_data)[1]
    #print(np.shape(band_data))
    add_para = 0
    print('E_shift is: ', E_shift)

    for k in range(4, rows_num_data):
        k += add_para
        if (band_data[0, k] == band_data[1, k] == 0 + E_shift) or \
                (band_data[0, k] == band_data[1, k] == 2 + E_shift):
            band_data = np.delete(band_data, k, axis=1)
            add_para -= 1
        else:
            print(band_data[0:5,k],file=f1)
    print(np.shape(band_data))
    f.close()

    columns_names = np.linspace(0, np.shape(band_data)[1] - 1, np.shape(band_data)[1])
    if os.path.exists(path + '/band_aims' + str(round(gap,2)) + '.dat'):
        pd_data = pd.read_csv(path + '/band_aims' + str(round(gap,2)) + '.dat')
        add_num *= 10
        columns_names += add_num
        start = np.shape(pd_data)[0]
        for i in range(np.shape(band_data)[0]):
            #print(np.shape(pd_data),np.shape(band_data))
            #print(start, i)
            pd_data.loc[start + i] = band_data[i]
        pd_data.to_csv(path + '/band_aims' + str(round(gap,2)) + '.dat', index=False)
    else:
        print('new')
        pd_data = pd.DataFrame(band_data, columns=columns_names.tolist())
        pd_data.to_csv(path + '/band_aims' + str(round(gap,2)) + '.dat', index=False)
pd_data = pd_data.drop(columns=['0.0', '1.0', '2.0','3.0'],axis=1)
pd_data.to_csv(path + '/band_aims' + str(round(gap,2)) + '.dat')

#-------------------------------get DOS data----------------------------------
f = open(path + '/KS_DOS_total_tetrahedron.dat')
lines = f.readlines()
raws_num = len(lines)  # 行数
cols_num = 2  # 列数
dos_data = np.zeros([raws_num-3, cols_num])
for i in range(raws_num-3):
    line = re.findall(r'.[0-9]*[.0-9][0-9]*', lines[i + 3])  # .[0-9]*[.0-9][0-9]*
    if line != []:
        dos_data[i, 0] = float(line[0]) + float(E_shift)
        dos_data[i, 1] = float(line[1])
    else:
        print('All words')
dos = pd.DataFrame(dos_data,index=None)
dos.to_csv(path + '/dos.dat')

plt.figure(dpi=250)
grid = plt.GridSpec(1, 3, wspace=0.1, hspace=0.1)

# ----------------------------Plot band------------------------------
fs = 10
lw = 0.5
band_file = '/band_aims' + str(round(gap,2)) + '.dat'
lines = pd.read_csv(path + band_file)

if path_type == 'm':
    xtick = [0,31,62,93]
    xtick_label = ['G','K','M','G']
    nrow = 93
elif path_type == 's':
    xtick = [0,31,62,93,124,155,186,217]
    xtick_label = ['G','K','M','G','A','H','L','A']
    nrow = 217
elif path_type == 'i':
    xtick = [0,31,62,93,124,155]
    xtick_label = ['G','K','M','G|R','G','U']
    nrow = 155

ncol = np.shape(lines)[1]
plt.subplot(grid[0, :2])
print('nrow is: ', nrow)
plt.axhline(y=0, xmin=0, xmax=nrow, ls=":",color='r',linewidth=lw)
for i in (xtick):
    plt.axvline(x=i, ls=":",color='k',linewidth=lw)

x = lines[lines.columns[0]][0:93]
for i in range(1,ncol):
    y = lines[lines.columns[i]][0:93]
    plt.plot(x,y, color='b',linewidth=0.7)
    #print(y)

if path_type == 's':
    x = lines[lines.columns[0]].loc[92:216]
    for i in range(1,ncol):
        y = lines[lines.columns[i]].loc[92:216]
        plt.plot(x,y, color='b',linewidth=0.7)
elif path_type == 'i':
    x = lines[lines.columns[0]].loc[93:154]
    for i in range(1,ncol):
        y = lines[lines.columns[i]].loc[93:154]
        plt.plot(x,y, color='b',linewidth=0.7)

plt.axis([0,nrow,E_min,E_max])
plt.ylabel("E - ${E_F}$ [eV]", fontsize=fs)
plt.xticks(xtick, xtick_label,fontsize=fs)
plt.yticks(np.arange(E_min,E_max,2),fontsize=fs)

# ----------------------------Plot DOS------------------------------
plt.subplot(grid[0,2])
plt.axhline(y=0, xmin=0, xmax=93, ls=":",color='r',linewidth=lw)
dos_filename = 'dos.dat'
f_dos = pd.read_csv(path + '/dos.dat')
x = f_dos['1']
y = f_dos['0']
fig = plt.plot(x,y,color='b',linewidth=0.7)
nrow = max(x)
plt.axis([0,nrow,E_min,E_max])
xtick_dos = set_xtick_dos(nrow)
plt.xticks(np.arange(0,nrow + 10,xtick_dos),fontsize=fs)
plt.yticks(np.arange(E_min,E_max,2), fontsize=fs)
ax = plt.gca()
ax.axes.yaxis.set_ticklabels([])

plt.savefig(path + '/band_dos_' + str(round(gap,2)) + '.png')
#plt.show()
f1.close()
f2.close()