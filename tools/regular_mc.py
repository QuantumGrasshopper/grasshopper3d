import numpy as np
from analytics import exact_ball_prob, optimal_r_inner, r0, r_outer, shell_prob

# set explicit shape: either ball or shell
# perform regular MC to get probability

def main():
    totalsteps = 100000
    counter = 0
    success = 0
    rng = np.random.default_rng(12345)

    jump = 1.2*r0
    r_inner = optimal_r_inner(jump).x[0]
    #r_inner = max(jump - r0, 0)
    outer_radius = r_outer(r_inner)

    for _ in range(totalsteps):

        # pick first point uniformly over volume using rejection method
        x = rng.uniform(-outer_radius,outer_radius)
        y = rng.uniform(-outer_radius,outer_radius)
        z = rng.uniform(-outer_radius,outer_radius)
        r = np.sqrt(x*x+y*y+z*z)

        if (r <= outer_radius) and (r >= r_inner):

            counter += 1

            # pick uniform point on sphere
            theta = rng.uniform(0,2*np.pi)
            u = rng.uniform(-1,1)
            relx = jump*np.sqrt(1-u*u)*np.cos(theta)
            rely = jump*np.sqrt(1-u*u)*np.sin(theta)
            relz = jump*u

            # displace grasshopper
            x = x+relx
            y = y+rely
            z = z+relz
            r = np.sqrt(x*x+y*y+z*z)

            # check success
            if (r <= outer_radius) and (r >= r_inner):
                success += 1

    print("Jump = ", jump)
    print("Inner radius = ", r_inner)
    print("The MC success probability equals: ", success/counter)
    print("The exact solid ball probability equals: ", exact_ball_prob(jump))
    print("The exact shell probability equals: ", shell_prob(jump, r_inner))


if __name__ == "__main__":
    main()
