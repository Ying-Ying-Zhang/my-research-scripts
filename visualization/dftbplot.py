#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2021/01/11
# @Function: plot band structure from dftb+ output

import numpy as np
import matplotlib.pyplot as plt
import os

#os.system(dp_bands band.out band)
#os.system(dp_dos band.out dos_tot.dat)

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

filename = 'band_tot.dat'
E_f = os.popen("grep 'Fermi level' detailed.out | awk '{print $5}'")
Ef = (E_f.read()).strip()
fs = 10
lw = 0.5
xtick_dos = 20
path_type = 'm' # monolayer; staggered; incliend

if path_type == 'm':
    xtick = [1,21,41,61]
    xtick_label = ['G','K','M','G']
    nrow = 61
elif path_type == 's':
    xtick = [1,21,41,61,101,121,141,161]
    xtick_label = ['G','K','M','G','A','H','L','A']
    nrow = 161
elif path_type == 'i':
    xtick = [1,21,41,61,101,141]
    xtick_label = ['G','K','M','G|R','G','U']
    nrow = 141

plt.figure(dpi=250)
grid = plt.GridSpec(1, 3, wspace=0.1, hspace=0.1)

f = np.loadtxt(filename)
ncol = np.shape(f)[1]
plt.subplot(grid[0, :2])
plt.axhline(y=float(Ef), xmin=0, xmax=nrow, ls=":",color='r',linewidth=lw)
for i in (xtick):
    plt.axvline(x=i, ls=":",color='k',linewidth=lw)
for i in range(1,ncol):
    x = f[0:61, 0]
    y = f[0:61,i]
    plt.plot(x,y, color='b',linewidth=0.7)
#print(np.shape(f),ncol)
if path_type == 's' or 'i':
    for i in range(1,ncol):
        x = f[61:161,0]
        y = f[61:161,i]
        plt.plot(x,y, color='b',linewidth=0.7)
plt.axis([1,nrow,-8,-2])
plt.ylabel("E [eV]", fontsize=fs)
plt.xticks(xtick, xtick_label,fontsize=fs)
plt.yticks(np.arange(-8,-2,2),fontsize=fs)

plt.subplot(grid[0,2])
plt.axhline(y=Ef, xmin=0, xmax=161, ls=":",color='r',linewidth=lw)
dos_filename = 'dos_tot.dat'
x = np.loadtxt(dos_filename, usecols=(1))
y = np.loadtxt(dos_filename, usecols=(0))
fig = plt.plot(x,y,color='b',linewidth=0.7)
plt.axis([0,max(x),-8,-2])
nrow = max(x)
xtick_dos = set_xtick_dos(nrow)
print(nrow, xtick_dos)
plt.xticks(np.arange(0,max(x),xtick_dos),fontsize=fs)
plt.yticks(np.arange(-8,-2,2), fontsize=fs)
ax = plt.gca()
ax.axes.yaxis.set_ticklabels([])

#plt.show()
plt.savefig('band_dos.png')
