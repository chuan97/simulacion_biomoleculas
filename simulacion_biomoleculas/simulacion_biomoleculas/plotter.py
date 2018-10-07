#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 13:08:24 2018

@author: juan
"""

import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')

f = open("/Users/juan/Xcode_projects/simulacion_biomoleculas/simulacion_biomoleculas/simulacion_biomoleculas/trajectoryVerlet.out", 'r')

t = [];
x = [];

for row in f:
    x.append(float(row.split(sep = ' ')[0]))
    t.append(float(row.split(sep = ' ')[1]))

plt.plot(x, t, "r")
plt.show()


