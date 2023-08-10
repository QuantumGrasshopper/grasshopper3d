import numpy as np
import matplotlib.pyplot as plt
from analytics import *

# plot optimal shell radius and corresponding probability
# compare with values for inner radius = jump - r0

jump_range = np.linspace(r0,2*r0,50)
opt_radii = []
opt_probs = []
reg_shell_probs = []
for i in jump_range:
    sol = optimal_r_inner(i)
    opt_radii.append(sol.x)
    opt_probs.append(-sol.fun)
    reg_shell_probs.append(shell_prob(i, i-r0))
    
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


