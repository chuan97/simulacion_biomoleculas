#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 13:08:24 2018

@author: juan
"""

import sys
import matplotlib.pyplot as plt
import itertools
import seaborn as sns
plt.style.use('seaborn-paper')

palette = itertools.cycle(sns.color_palette())

filename_1 = sys.argv[1]
filename_2 = sys.argv[2]
filename_3 = sys.argv[3]
filename_4 = sys.argv[4]

f_1 = open(filename_1, 'r')
f_2 = open(filename_2, 'r')
f_3 = open(filename_3, 'r')
f_4 = open(filename_4, 'r')

t = [];
x_1 = [];
x_2 = [];
x_3 = [];
x_4 = [];

for row in f_1:
    t.append(float(row.split(sep = ' ')[0]))
    x_1.append(float(row.split(sep = ' ')[1]))

for row in f_2:
    x_2.append(float(row.split(sep = ' ')[1]))

for row in f_3:
    x_3.append(float(row.split(sep = ' ')[1]))

for row in f_4:
    x_4.append(float(row.split(sep = ' ')[1]))

plt.subplot(2, 2, 1)
plt.plot(t, x_1, color = next(palette))
plt.title('damping = 0.0')
plt.ylabel('x')

plt.subplot(2, 2, 2)
plt.plot(t, x_2, color = next(palette))
plt.title('damping = 0.1')

plt.subplot(2, 2, 3)
plt.plot(t, x_3, color = next(palette))
plt.title('damping = 1.0')
plt.xlabel('t')
plt.ylabel('x')

plt.subplot(2, 2, 4)
plt.plot(t, x_4, color = next(palette))
plt.title('damping = 10.0')
plt.xlabel('t')

plt.tight_layout()
plt.suptitle('Euler Maruyama')
plt.subplots_adjust(top=0.88)

plt.show()


