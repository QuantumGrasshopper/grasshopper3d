#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as anim

# generate animation of simulated annealing using saved configuration sequence

def read_result_value(label, value_type):
    with open("result.dat") as f:
        for line in f:
            key, separator, value = line.partition(":")
            if separator and key.strip() == label:
                return value_type(value.strip())
    raise ValueError(f"Could not find '{label}' in result.dat")

gridsize = read_result_value("Size of grid", int)

data = np.atleast_2d(np.loadtxt("config.dat", dtype=float))

fig = plt.figure()
ax = plt.axes(projection='3d',elev=30, azim=45)

numframes = len(data)

def animation_function(i):
    ax.cla()
    row = data[i]
    coordinates = row[:-1].astype(int)    # last element is raw energy
    z = coordinates//gridsize//gridsize
    y = coordinates//gridsize - z*gridsize
    x = coordinates - y*gridsize - z*gridsize*gridsize
    
    ax.set_xlim3d(0,gridsize)
    ax.set_ylim3d(0,gridsize)
    ax.set_zlim3d(0,gridsize)
    ax.set_box_aspect((1,1,1))
    ax.view_init(elev=30, azim=45)
    ax.scatter(x, y, z, c=z, cmap='viridis')

  
animation = anim.FuncAnimation(fig, animation_function, 
                               frames=numframes, interval = 2)

plt.show()
