#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d
import matplotlib.animation as anim
import plotly.graph_objects as go

# generate animation of simulated annealing using saved configuration sequence

data = np.genfromtxt("config.dat", dtype = int)

with open("result.dat") as f:
    lines = f.readlines()
    gridsize = int(lines[4].strip().split(" ")[-1])

fig = plt.figure()
ax = plt.axes(projection='3d',elev=30, azim=45)

numframes = len(data)

def animation_function(i):
    ax.cla()
    row = data[i,:-1]    # last element is energy
    z = row//gridsize//gridsize
    y = row//gridsize - z*gridsize
    x = row - y*gridsize - z*gridsize*gridsize
    
    ax.set_xlim3d(0,gridsize)
    ax.set_ylim3d(0,gridsize)
    ax.set_zlim3d(0,gridsize)
    ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
    ax.view_init(elev=30, azim=45)
    ax.scatter(x, y, z, c=z, cmap='viridis')

  
animation = anim.FuncAnimation(fig, animation_function, 
                               frames=numframes, interval = 2)

plt.show()
