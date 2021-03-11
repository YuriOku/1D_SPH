# 1D_SPH
This is a one-dimensional smoothed particle hydrodynamic simulation code written in Julia.
The Sod shock tube problem is generated and solved by the SPH scheme.

# Usage
1. edit main.jl to setup the problem
2. (optional) `export JULIA_NUM_THREAD=8`
3. `julia main.jl`

# References
## review papers on SPH scheme
- Rosswog, S., ``Astrophysical smooth particle hydrodynamics'', <i>New Astronomy Reviews</i>, vol. 53, no. 4–6, pp. 78–104, 2009. doi:10.1016/j.newar.2009.08.007.
- Rosswog, S., ``SPH Methods in the Modelling of Compact Objects'', <i>Living Reviews in Computational Astrophysics</i>, vol. 1, no. 1, 2015. doi:10.1007/lrca-2015-1.
- Springel, V., ``Smoothed Particle Hydrodynamics in Astrophysics'', <i>Annual Review of Astronomy and Astrophysics</i>, vol. 48, pp. 391–430, 2010. doi:10.1146/annurev-astro-081309-130914.

## SPH formulation
- Saitoh, T. R. and Makino, J., “A Density-independent Formulation of Smoothed Particle Hydrodynamics”, <i>The Astrophysical Journal</i>, vol. 768, no. 1, 2013. doi:10.1088/0004-637X/768/1/44.
- Rosswog, S., “The Lagrangian hydrodynamics code MAGMA2”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 498, no. 3, pp. 4230–4255, 2020. doi:10.1093/mnras/staa2591.
- Wadsley, J. W., Keller, B. W., and Quinn, T. R., “Gasoline2: a modern smoothed particle hydrodynamics code”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 471, no. 2, pp. 2357–2369, 2017. doi:10.1093/mnras/stx1643.

## kernel gradients
- García-Senz, D., Cabezón, R. M., and Escartín, J. A., “Improving smoothed particle hydrodynamics with an integral approach to calculating gradients”, <i>Astronomy and Astrophysics</i>, vol. 538, 2012. doi:10.1051/0004-6361/201117939.
- Rosswog, S., “Boosting the accuracy of SPH techniques: Newtonian and special-relativistic tests”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 448, no. 4, pp. 3628–3664, 2015. doi:10.1093/mnras/stv225.

## kernel functions
- Dehnen, W. and Aly, H., ``Improving convergence in smoothed particle hydrodynamics simulations without pairing instability'', <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 425, no. 2, pp. 1068–1082, 2012. doi:10.1111/j.1365-2966.2012.21439.x.

## time integration scheme
- Springel, V., Pakmor, R., Zier, O., and Reinecke, M., “Simulating cosmic structure formation with the GADGET-4 code”, <i>arXiv e-prints</i>, 2020.
- Ketcheson, D. I., “Runge-Kutta methods with minimum storage implementations”, <i>Journal of Computational Physics</i>, vol. 229, no. 5, pp. 1763–1773, 2010. doi:10.1016/j.jcp.2009.11.006.
- Stone, J. M., Tomida, K., White, C. J., and Felker, K. G., “The Athena++ Adaptive Mesh Refinement Framework: Design and Magnetohydrodynamic Solvers”, <i>The Astrophysical Journal Supplement Series</i>, vol. 249, no. 1, 2020. doi:10.3847/1538-4365/ab929b.

## exact Riemann solver
- Toro, E. F., ``Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction'', Springer, 2012. doi:10.1007/b79761
