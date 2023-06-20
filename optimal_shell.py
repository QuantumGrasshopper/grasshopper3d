import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.optimize as opt 

r0 = (3/4/np.pi)**(1/3)

def ball_with_hole_prob(jump):

    # this is for a specific r_inner, typically d-r0
    r_inner = max(jump - r0, 0)
    r_outer = (3/4/np.pi+r_inner**3)**(1/3)
    
    def integrand1(r):
        return 2*np.pi*(r*r+r*(r_outer**2-jump**2-r**2)/2/jump)
    def integrand2(r):
        return np.pi*r*(r_outer**2-r_inner**2)/jump
    
    return quad(integrand1, r_inner, jump-r_inner)[0]+quad(integrand2, jump-r_inner, r_outer)[0]

def optimal_r_inner(jump):

    def r_outer(r_inner):
        return (3/4/np.pi+r_inner**3)**(1/3)
    
    def integrand1(r, r_inner):
        return 2*np.pi*(r*r+r*(r_outer(r_inner)**2-jump**2-r**2)/2/jump)
    def integrand2(r, r_inner):
        return 2*np.pi*(r*r+r*(r_inner**2-jump**2-r**2)/2/jump)
    def integral(r_inner):
        return quad(integrand1, r_inner, r_outer(r_inner), args=(r_inner,))[0]-quad(integrand2, jump-r_inner, r_outer(r_inner), args=(r_inner,))[0]
    
    def negintegral(r_inner):
        return -integral(r_inner)
    
    return opt.minimize(negintegral, jump-r0)
    
jump_range = np.linspace(r0,2*r0,50)
opt_radii = []
opt_probs = []
reg_shell_probs = []
for i in jump_range:
    sol = optimal_r_inner(i)
    opt_radii.append(sol.x)
    opt_probs.append(-sol.fun)
    reg_shell_probs.append(ball_with_hole_prob(i))
    
plt.plot(jump_range,opt_radii,label=r'$R_i$ optimal')
plt.plot(jump_range,jump_range-r0, label=r'$R_i=d-r_0$')
plt.xlim(r0,2*r0)
plt.ylim(0,0.7)
plt.xlabel(r'$d$')
plt.ylabel('hole radius $R_i$')
plt.legend(loc='upper left')
plt.show()

plt.plot(jump_range,opt_probs,label=r'$R_i$ optimal')
plt.plot(jump_range,reg_shell_probs,label=r'$R_i=d-r_0$')
plt.xlim(r0,2*r0)
plt.ylim(0,0.35)
plt.xlabel(r'$d$')
plt.ylabel("probability")
plt.legend()
plt.show()


