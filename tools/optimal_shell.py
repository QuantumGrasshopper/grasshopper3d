# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Olga Goulko and David Llamas

import numpy as np
import matplotlib.pyplot as plt
from analytics import optimal_r_inner, r0, shell_prob

# plot optimal shell radius and corresponding probability
# compare with values for inner radius = jump - r0

def main():
    jump_range = np.linspace(r0,2*r0,50)
    opt_radii = []
    opt_probs = []
    reg_shell_probs = []
    for jump in jump_range:
        solution = optimal_r_inner(jump)
        opt_radii.append(solution.x[0])
        opt_probs.append(-solution.fun)
        reg_shell_probs.append(shell_prob(jump, jump-r0))

    plt.figure()
    plt.plot(jump_range,opt_radii,label=r'$R_i$ optimal')
    plt.plot(jump_range,jump_range-r0, label=r'$R_i=d-r_0$')
    plt.xlim(r0,2*r0)
    plt.ylim(0,0.7)
    plt.xlabel(r'$d$')
    plt.ylabel('hole radius $R_i$')
    plt.legend(loc='upper left')
    plt.show()

    plt.figure()
    plt.plot(jump_range,opt_probs,label=r'$R_i$ optimal')
    plt.plot(jump_range,reg_shell_probs,label=r'$R_i=d-r_0$')
    plt.xlim(r0,2*r0)
    plt.ylim(0,0.35)
    plt.xlabel(r'$d$')
    plt.ylabel("probability")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
