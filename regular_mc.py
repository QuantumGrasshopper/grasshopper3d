import numpy as np
import matplotlib.pyplot as plt

# Setting explicit shape and doing regular MC to get probability

totalsteps = 1000000
r0 = (3/4/np.pi)**(1/3)
jump = 0.3
counter = 0
success = 0

def exact_ball_prob(jump):
    return 1 - np.pi*r0*r0*jump + np.pi/12*jump*jump*jump

for steps in range(totalsteps):

    # shape is a ball
    # pick first point uniformly over volume using rejection method

    x = np.random.uniform(-r0,r0);
    y = np.random.uniform(-r0,r0);
    z = np.random.uniform(-r0,r0);
    r = np.sqrt(x*x+y*y+z*z)
    
    if r <= r0:
        
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
        if r <= r0:
            success += 1
 
print("The MC success probability equals: ", success/counter)
print("The exact probability equals: ", exact_ball_prob(jump))
        
        
