import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize

r0 = (3/4/np.pi)**(1/3)

def exact_ball_prob(jump):
    """Return the unit-volume ball probability for 0 <= jump <= 2*r0."""
    jump = float(jump)
    if not np.isfinite(jump) or jump < 0 or jump > 2*r0:
        raise ValueError(f"Jump must be finite and in [0, {2*r0}].")
    return 1 - np.pi*r0*r0*jump + np.pi/12*jump*jump*jump

def r_outer(r_inner):
    """Return the outer radius of a unit-volume shell with inner radius r_inner."""
    if not np.all(np.isfinite(r_inner)) or np.any(np.asarray(r_inner) < 0):
        raise ValueError("Inner radius must be finite and non-negative.")
    return (3/4/np.pi+r_inner**3)**(1/3)

def shell_prob(jump, r_inner):
    
    def integrand1(r):
        return 2*np.pi*(r*r+r*(r_outer(r_inner)**2-jump**2-r**2)/2/jump)
    def integrand2(r):
        return np.pi*r*(r_outer(r_inner)**2-r_inner**2)/jump
    
    return quad(integrand1, r_inner, jump-r_inner)[0] +\
           quad(integrand2, jump-r_inner, r_outer(r_inner))[0]

def optimal_r_inner(jump):
   
    def integrand1(r, r_inner):
        return 2*np.pi*(r*r+r*(r_outer(r_inner)**2-jump**2-r**2)/2/jump)
    def integrand2(r, r_inner):
        return 2*np.pi*(r*r+r*(r_inner**2-jump**2-r**2)/2/jump)
    def integral(r_inner):
        return quad(integrand1, r_inner, r_outer(r_inner), args=(r_inner,))[0] -\
               quad(integrand2, jump-r_inner, r_outer(r_inner), args=(r_inner,))[0]
    
    def negintegral(r_inner):
        return -integral(r_inner)
    
    return minimize(negintegral, jump-r0)
