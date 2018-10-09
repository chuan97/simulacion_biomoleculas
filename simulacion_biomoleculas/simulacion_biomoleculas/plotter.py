#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 13:08:24 2018

@author: juan
"""

import sys
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')

filename = sys.argv[1]

f = open(filename, 'r')

t = [];
x = [];
v = [];
V = [];
T = [];
E = [];

for row in f:
    t.append(float(row.split(sep = ' ')[0]))
    x.append(float(row.split(sep = ' ')[1]))
    v.append(float(row.split(sep = ' ')[2]))
    V.append(float(row.split(sep = ' ')[3]))
    T.append(float(row.split(sep = ' ')[4]))
    E.append(float(row.split(sep = ' ')[5]))
    
    

plt.plot(t, V)
plt.plot(t, T)
plt.plot(t, E)
plt.show()


