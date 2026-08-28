#!/usr/bin/python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Olga Goulko

import numpy as np
import matplotlib.pyplot as plt

# quick plot saved spin configuration
# no transparency, cannot see inside configuration

def read_result_value(label, value_type):
    with open("result.dat") as f:
        for line in f:
            key, separator, value = line.partition(":")
            if separator and key.strip() == label:
                return value_type(value.strip())
    raise ValueError(f"Could not find '{label}' in result.dat")

gridsize = read_result_value("Size of grid", int)

inputfile = input("Name of file containing spin configuration: ")
data = np.loadtxt(inputfile, dtype = int)

z = data//gridsize//gridsize
y = data//gridsize - z*gridsize
x = data - y*gridsize - z*gridsize*gridsize


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim3d(0,gridsize)
ax.set_ylim3d(0,gridsize)
ax.set_zlim3d(0,gridsize)
ax.set_box_aspect((1,1,1))
ax.scatter(x, y, z, c=z, cmap='viridis')

plt.show()
