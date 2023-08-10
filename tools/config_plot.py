#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# quick plot saved spin configuration
# no transparency, cannot see inside configuration

inputfile = input("Name of file containing spin configuration: ")

data = np.genfromtxt(inputfile, dtype = int)

with open("result.dat") as f:
    lines = f.readlines()
    gridsize = int(lines[4].strip().split(" ")[-1])

z = data//gridsize//gridsize
y = data//gridsize - z*gridsize
x = data - y*gridsize - z*gridsize*gridsize


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x, y, z, c=z, cmap='viridis')

plt.show()
