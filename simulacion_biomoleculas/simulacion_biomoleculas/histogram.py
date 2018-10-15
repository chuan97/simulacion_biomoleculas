#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 09:30:08 2018

@author: juan
"""

import sys
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_palette(sns.color_palette("hls", 20))

filename = sys.argv[1]
f = open(filename, 'r')

x = []
y = []

for row in f:
    x.append(float(row.split(sep = ' ')[0]))
    y.append(float(row.split(sep = ' ')[1]))
    
plt.plot(x, y)
plt.show()
