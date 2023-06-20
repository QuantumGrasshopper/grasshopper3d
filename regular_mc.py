import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.optimize as opt

# Setting explicit shape and doing regular MC to get probability

totalsteps = 1000000
counter = 0
success = 0

r0 = (3/4/np.pi)**(1/3)

def exact_ball_prob(jump):
    return 1 - np.pi*r0*r0*jump + np.pi/12*jump*jump*jump

def ball_with_hole_prob(jump):

    r_inner = max(jump - r0, 0)
    r_outer = (3/4/np.pi+r_inner**3)**(1/3)
    
    def integrand1(r):
        return 2*np.pi*(r*r+r*(r_outer**2-jump**2-r**2)/2/jump)
    def integrand2(r):
        return np.pi*r*(r_outer**2-r_inner**2)/jump
    
    return quad(integrand1, r_inner, jump-r_inner)[0]+quad(integrand2, jump-r_inner, r_outer)[0]
    # return quad(integrand1, r_inner, r0)[0]+quad(integrand2, r0, r_outer)[0]
  

    
jump = 1.2*r0
r_inner = max(jump - r0, 0)
r_outer = (3/4/np.pi+r_inner**3)**(1/3)

for steps in range(totalsteps):

    # shape is a ball (possibly with hole)
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
 
print("The MC success probability equals: ", success/counter)
print("The exact solid ball probability equals: ", exact_ball_prob(jump))
print("The exact shell probability equals: ", ball_with_hole_prob(jump))
        
