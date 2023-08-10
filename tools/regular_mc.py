import numpy as np
from analytics import *

# set explicit shape: either ball or shell
# perform regular MC to get probability

totalsteps = 100000
counter = 0
success = 0

jump = 1.2*r0
r_inner = optimal_r_inner(jump).x[0]
#r_inner = max(jump - r0, 0)
r_outer = r_outer(r_inner)

for steps in range(totalsteps):

    # pick first point uniformly over volume using rejection method
    x = np.random.uniform(-r_outer,r_outer);
    y = np.random.uniform(-r_outer,r_outer);
    z = np.random.uniform(-r_outer,r_outer);
    r = np.sqrt(x*x+y*y+z*z)
    
    if (r <= r_outer) and (r >= r_inner):
        
        counter += 1
        
        # pick uniform point on sphere
        theta = np.random.uniform(0,2*np.pi)
        u = np.random.uniform(-1,1)
        relx = jump*np.sqrt(1-u*u)*np.cos(theta)
        rely = jump*np.sqrt(1-u*u)*np.sin(theta)
        relz = jump*u
        
        # displace grasshopper
        x = x+relx
        y = y+rely
        z = z+relz
        r = np.sqrt(x*x+y*y+z*z)
        
        # check success
        if (r <= r_outer) and (r >= r_inner):
            success += 1


print("Jump = ", jump)
print("Inner radius = ", r_inner)
print("The MC success probability equals: ", success/counter)
print("The exact solid ball probability equals: ", exact_ball_prob(jump))
print("The exact shell probability equals: ", shell_prob(jump, r_inner))
