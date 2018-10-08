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

f = open("/Users/juan/Xcode_projects/simulacion_biomoleculas/simulacion_biomoleculas/DerivedData/simulacion_biomoleculas/Build/Products/Debug/" + filename + ".out", 'r')

t = [];
x = [];

for row in f:
    x.append(float(row.split(sep = ' ')[1]))
    t.append(float(row.split(sep = ' ')[0]))

plt.plot(t, x)
plt.show()


