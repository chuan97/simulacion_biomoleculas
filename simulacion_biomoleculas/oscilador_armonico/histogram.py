#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 09:30:08 2018

@author: juan
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_palette(sns.color_palette("hls", 20))

filenames = [sys.argv[1],
             sys.argv[2],
             sys.argv[3],
             sys.argv[4],
             sys.argv[5],
             sys.argv[6]]

for name, i in zip(filenames, range(6)):
    x = []
    y = []
    
    f = open(name, 'r')
    
    for row in f:
        x.append(float(row.split(sep = ' ')[0]))
        y.append(float(row.split(sep = ' ')[1]))

    if i < 3:
        plt.subplot(1, 2, 1)
        end = name.find('.out')
        beg = name.find('_', end - 5) + 1
        plt.plot(x, y, label = name[beg : end])
    

    else:
        plt.subplot(1, 2, 2)
        end = name.find('.out')
        beg = name.find('_', end - 5) + 1
        plt.plot(x, y, label = name[beg : end])

plt.subplot(1, 2, 1)
plt.plot(x, [1 / np.sqrt(2 * np.pi) * np.exp(- 0.5 * t ** 2) for t in x], label = 'guía', zorder = 0)
plt.legend(title = 'damping', loc = 'upper right', frameon = False)
plt.title('Posiciones')
plt.xlabel('x')

plt.subplot(1, 2, 2)
plt.plot(x, [1 / np.sqrt(2 * np.pi) * np.exp(- 0.5 * t ** 2) for t in x], label = 'guía', zorder = 0)
plt.legend(title = 'damping', loc = 'upper right', frameon = False)
plt.title('Velocidades')
plt.xlabel('v')

plt.tight_layout()
first_ = filenames[0].find('_')
second_ = filenames[0].find('_', first_ + 1)
third_ =filenames[0].find('_', second_ + 1)
plt.suptitle(filenames[0][second_ + 1 : third_].upper())
plt.subplots_adjust(top=0.88)
plt.show()

