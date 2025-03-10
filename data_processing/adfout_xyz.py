#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 2021/02/22
# @Function: extract optimized structure from adf output

import re
import os
import fnmatch


path = os.getcwd()

pattern_1 = re.compile(r'.*Atoms')  # end of head
pattern_2 = re.compile(r'End')  # start of tail
pattern_3 = re.compile(r'.*Index Symbol')  # start of coordination
pattern_4 = re.compile(r'Lattice vectors (angstrom)*')
pattern_5 = re.compile('.*BondOrders')

line_num_Atoms = 0
line_num_End = 0
line_num_IndexSym = 0
file_run = fnmatch.filter(os.listdir(path),'*.run')[0]
file_out = fnmatch.filter(os.listdir(path),'*.out')[0]
print(file_run,file_out)
with open(file_run,"r") as frun, open(file_out,"r") as fout:
    freadrun = frun.readlines()
    freadout = fout.readlines()
    for cnt, line in enumerate(freadrun):
        if re.match(pattern_1,line):  # get the index of the end of head
            line_num_Atoms=cnt+1
        if re.match(pattern_5,line):
            line_num_Bondsorder=cnt
        elif re.match(pattern_2,line):     # get the index of the start of tail
            line_num_End = cnt+1
            break
    for cnt, line in enumerate(freadout):   # get the index of the start of coordination in .out file
        if re.match(pattern_3,line):
            line_num_IndexSym=cnt+1
        if re.match(pattern_4, line):
            line_num_lattice=cnt+1

    with open("new_run","w") as fw, open ("geo_end.xyz","w") as geo:
        atom_num = 0
        geo.write('\n\n\n')
        for line in freadrun[0:line_num_Atoms]:
            fw.write(str(line))
        for line in freadout[line_num_IndexSym:]:
            if re.match(r'\s*\n',line):
                break
            else:
                atom_num += 1
                fw.write(str(re.sub(r'\s*\d*',' ',line, 1)))
                geo.write(str(re.sub(r'\s*\d*',' ',line, 1)))
        
        if line_num_Atoms > 7:  # adf/2017
            fw.write('    End\n\n    Charge 0\n\n    Lattice\n')
            for line in freadout[line_num_lattice:line_num_lattice + 3]:
                if re.match(r'\s*\d\s*\d*.+\s*', line):
                    fw.write(str(re.sub(r'\s*\d*', ' ', line, 1)))
                    vec = 'VEC' + str(re.split(r'\s*', line)[2])
                    geo.write(str(re.sub(r'\s*\d*', vec, line, 1)))
            for line in freadrun[line_num_End - 2:]:
                fw.write(str(line))
        else:  # adf/2019
            fw.write('    End\n    Lattice\n')
            for line in freadout[line_num_lattice:line_num_lattice + 3]:
                if re.match(r'\s*\d\s*\d*.+\s*', line):
                    fw.write(str(re.sub(r'\s*\d*', ' ', line, 1)))
                    vec = 'VEC' + str(re.split(r'\s*', line)[2])
                    geo.write(str(re.sub(r'\s*\d*', vec, line, 1)))
            fw.write('    End\n')
            for line in freadrun[line_num_Bondsorder:]:
                fw.write(str(line))
        geo.seek(0, 0)
        geo.write(str(atom_num)+'\n\n')


