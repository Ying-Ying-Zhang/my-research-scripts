#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2021/03/11
# @Function: Simulate XRD and plot the figure

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import fnmatch

# calculate XRD, which would give us discrete data
def cal_xrd(file_name, two_theta_range=(0, 50)):
    fname, ext = os.path.splitext(file_name)
    print(file_name,type(file_name))
    structure = Structure.from_file(file_name)
    c = XRDCalculator(wavelength='CuKa1')
    xrd = XRDCalculator.get_pattern(c, structure, scaled=True, two_theta_range=two_theta_range)
    return xrd


# Convert discrete data into continuous data with equally spaced data by using Gaussian and Lorentzian fitting
def spectrum(x, y, sigma, x_range):
    gE = []
    for xi in x_range:
        tot = 0
        for xj, oss in zip(x, y):
            L = (FWHM / (2 * np.pi)) * (1 / ((xj - xi) ** 2 + 0.25 * FWHM ** 2))
            G = a * np.exp(-(xj - xi) ** 2 / (2 * sigma ** 2))
            P = omega * G + (1 - omega) * L
            tot += oss * P
            # tot+=oss*np.exp(-((((xj-xi)/sigma)**2)))
        gE.append(tot)
    return gE


# Plot the XRD figure
def plot_xrd(x, y, two_theta_range=(0, 50)):
    fig, ax = plt.subplots(figsize=(6, 4), dpi=250)
    ax.plot(x, y)
    ax.set_xlabel("2θ(degrees)", fontsize=16)
    ax.xaxis.set_tick_params(labelsize=14, width=1.5)
    ax.yaxis.set_tick_params(labelsize=14, width=1.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.set_xlim(two_theta_range)
    ax.set_ylabel("Intensity", fontsize=16)
    plt.tight_layout()
    plt.savefig(path + '/generate_xrd.png')

#path = '/home/yizh781b/COF-1/PXRD/mesitylene/'
path = '/PATH of the structures to be calculated/'
x_min, x_max = 0, 50  # the range of two_theta

# parameter settings for the fitting
FWHM = 0.1
sigma = FWHM * 0.42463
omega = 0.001  # the proportion of Gaussian fitting
a = 1/(sigma*np.sqrt(2*np.pi))
x_num = int((x_max - x_min) / 0.02) + 1
xrd_x = np.linspace(x_min, x_max, num=x_num, endpoint=True)
# list all the files with '.cif' in the specified directory, one can also change the format
files = fnmatch.filter(os.listdir(path),'*.cif')
os.chdir(path)

for file in files:
    print(file,type(file))
    f = pd.DataFrame(xrd_x, columns=['twotheta'])
    f.to_csv('pxrd.csv', index=False)
    fname, ext = os.path.splitext(file)
    print('calculating', str(file))
    xrd_data = cal_xrd(str(file), two_theta_range=(x_min, x_max))
    gxrd = spectrum(xrd_data.x, xrd_data.y, sigma, xrd_x)
    f['intensity_' + str(fname)] = gxrd
f.to_csv('pxrd.csv', index=False)  # write the XRD data to csv file


