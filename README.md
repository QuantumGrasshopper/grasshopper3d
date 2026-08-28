# Grasshopper3D

Grasshopper3D is a C++ implementation of simulated annealing for the three dimensional Euclidean grasshopper problem on a cubic grid. The lawn has unit volume and is represented by `N` occupied grid cells. The code searches for configurations that maximize the discretized grasshopper success probability for a fixed jump distance `d`.

## Citation

If you use Grasshopper3D in published work, please cite the software release (see [`CITATION.cff`](CITATION.cff)) and the following papers:

- D. Llamas, J. Kent-Dobias, K. Chen, A. Kent, and O. Goulko,
  *Origin of Symmetry Breaking in the Grasshopper Model*,
  Phys. Rev. Research **6**, 023235 (2024),
  https://doi.org/10.1103/PhysRevResearch.6.023235.

- O. Goulko and A. Kent, *The grasshopper problem*, 
  Proc. R. Soc. A **473**, 20170494 (2017),
  https://doi.org/10.1098/rspa.2017.0494.

## Funding

This work was supported by the National Science Foundation under Grant Nos. PHY-2112738 ("CQIS: The Grasshopper Problem") and OSI-2328774 ("ExpandQISE: Track 2: EQUIP-UMB").

## Related software

- The analogous implementation for the 2D grasshopper problem is available at https://github.com/QuantumGrasshopper/grasshopper2d
- A broader Python implementation supporting both spherical and Euclidean grasshopper models is available at https://github.com/llamas7/grasshopper.

## Prerequisites and compilation

- make
- GSL library
- C++17 compiler (GCC and Clang tested)

To compile the code, type `make` in the shell while in the folder containing the code. This uses the default GCC compiler. To use a different compiler, for example Clang, type `make CXX=clang++`, etc.

## Testing

Run the unit tests with

> `make test`

run the integration tests with

> `make integration-test`

or run the complete test suite with

> `make check`

## Running the code

The code uses the following command-line options:

| option          | comment |
| --------------- | ------------------ |
| `-d`            | required positive grasshopper jump distance in regular length units (unit in which the lawn volume is 1), e.g. 0.1 or 0.5 or 2 or 10 |
| `-N`            | total number of grid points with spin 1 (10000 is default) |
| `-gridsize`     | cubic-grid edge length in cells (automatically sized if omitted); must satisfy `floor((gridsize - 1)/2) >= ceil(d/cellSize) + 1`, where `cellSize=N^(-1/3)`; undersized grids are rejected |
| `-hours`        | how many hours the code should run (can be less than 1 hour), e.g. 0.1, 0.5, 2, 48, can also be 0 (default value) if you only want to look at the initial configuration |
| `-steps`        | how many steps the code should maximally run (`1e12` by default), but code will terminate earlier if maximal time is reached |
| `-tempsteps`    | initial number of steps before first temperature decrease (defaults to `N`; the number of steps between decreases goes up with each round) |
| `-inittemp`     | initial temperature (20 by default) |
| `-fintemp`      | final temperature (0.1 by default; need to run long enough to reach it) |
| `-annealsteps`  | number of simulated annealing steps between initial and final temperature (1000 by default) |
| `-initconf`     | how to initialise the system: currently implemented: `random` (default), or `load` (load configuration from file called `initconf.dat`) |
| `-delta`        | choice of delta-function discretization: exactly `0` (default) or `1`; see Goulko and Kent (2017) for their definitions |
| `-randomseed`   | unsigned initial value for the random number generator (if omitted or set to `0`, a seed is generated from the system clock) |
| `-overwrite`    | output overwrite policy: exactly `0` (default, reject if an output artifact exists) or `1` (remove old output artifacts before starting) |

- Required options: `d`
- Recommended options: `N`, `gridsize`, `hours`

Every option accepts exactly one value and may be supplied at most once. Unknown options, missing values, malformed or out-of-range numbers, and non-finite floating-point values are rejected. The number of spins must satisfy `0 < N < gridsize^3` after automatic or explicit grid sizing.

The code uses the following standard output files:

file name          | comment
---------------    | ---------
`result.dat`       | general info about the simulation and parameters
`initconf.dat`     | initial spin configuration; with `-initconf load`, this is an input and is never removed or overwritten
`finconf.dat`      | final spin configuration
`bestconf.dat`     | best configuration over the whole run
`energies.dat`     | every annealing round prints the current raw pair energy
`temperatures.dat` | every annealing round prints the counter, the current temperature, and the current acceptance ratio
`config.dat`       | stores selected spin configurations and raw energies across the annealing trajectory for animation or analysis; output is capped at 100 MB

By default, a run is rejected before creating files if any standard output file listed above already exists. With `-overwrite 1`, the complete standard output set is removed before the run. The only exception is `initconf.dat` with `-initconf load`, which is preserved as the run's input. A cleanup error aborts the run before simulation output begins.

Example of command to run the code:

> `./grasshopper -N 10000 -initconf random -gridsize 50 -d 0.3 -hours 1 -inittemp 20.0 -fintemp 0.05 `

## Python tools

Python scripts for plotting spin configurations and other tools can be found in the subfolder `tools`.

- `config_plot.py` plots a grasshopper spin configuration in 3D, such as `initconf.dat` or `bestconf.dat` etc. Input is the name of the configuration file; everything else is read automatically from `result.dat`, which must be present.
- `anneal_animation.py` generates a 3D animation of the simulated annealing process from the file `config.dat`, which is output by the main code. The size of the grid is read automatically from `result.dat`, which must be present.

Additional scripts provide analytical and Monte Carlo calculations used for comparisons with three-dimensional grasshopper configurations, see also Llamas et al. (2024):

- `analytics.py` implements analytical success probabilities for unit-volume balls and spherical shells, including optimization of the shell inner radius.
- `optimal_shell.py` plots the optimal shell inner radius and corresponding success probability as functions of jump distance, including comparison with the shell defined by `R_i=d-r_0`.
- `regular_mc.py` performs a direct Monte Carlo calculation for a representative optimized spherical shell and compares the result with the analytical ball and shell probabilities.

The visualization tools require NumPy and Matplotlib. The analytical tools additionally require SciPy.

## License

This software is distributed under the GNU General Public License version 3 or, at your option, any later version (GPL-3.0-or-later). See [LICENSE](LICENSE) for details.

Third-party components retain their respective licenses.
