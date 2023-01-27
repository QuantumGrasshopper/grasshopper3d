import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

inputfile = input("Name of file containing spin configuration: ")
gridsize = int(input("Grid size: "))

data = np.genfromtxt(inputfile, dtype = int)

numberspins = len(data)

z = data//gridsize//gridsize
y = data//gridsize - z*gridsize
x = data - y*gridsize - z*gridsize*gridsize


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x, y, z, c=z, cmap='viridis')

plt.show()
