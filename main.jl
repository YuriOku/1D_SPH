# simulation settings

# kernel function: "cubic spline", "Wendland C2"
kernel = "cubic spline"
# kernel = "Wendland C2"

# kernel interpolation: "standard"
gradient = "standard"

# sph formulation: "vanilla ice"
formulation = "vanilla ice"

# volume element: "mass", "U"
volume_element = "U"

# time integrator: "RK2", "VL2", "leap frog", "symplectic Euler"
time_integrator = "leap frog"

# output file name
outputfile = "out/out"

# The Courant-Friedrichs-Levi (CFL) number
cfl = 0.3

# time-dependent viscosity
time_dependent_viscosity = false
alpha_max = 1.5
alpha_min = 0.1

# upper limit of time step
dt_max = 1

# factor for evaluation of smoothing length
# this code uses Gadget's definition of smoothing length (smoothing kernel vanishes at 1h rather than at 2h)
eta_hsml = 2.4

# maximum number of SPH neighbour particles
Nngb = 32

# simulation time and output timing
t_start = 0
t_end = 0.2
output_interval = 0.001

# specific heat ratio
gamma = 1.4

# number of particles
Npart = 1000

# simulation box boundaries
x_min = -1
x_max = 1
center = 0

# density and pressure in the left side
rho_left = 1
v_left = 0
P_left = 1

# density and pressure in the right side
rho_right = 0.25
v_right = 0
P_right = 0.1795

# plot density, pressure and velocity: true, false
plot_figure = true
x_min_plot = -0.5
x_max_plot = 0.5

# number of sample points for exact Riemann solver 
Nsample_riemann = 1e3
# tolerance criterion in iterative process in exact Riemann solver
TOL = 1e-6

# check conservation of momentum and energy
debug = true

# calculation

if plot_figure
  using GR
  using Plots
  ENV["GKSwstype"] = "100"
  gr()
end
using Printf

include("./src/begin.jl")
include("./src/evaluate.jl")
include("./src/exact_riemann.jl")
include("./src/finish.jl")
include("./src/force.jl")
include("./src/integrate.jl")
include("./src/kernel.jl")
include("./src/output.jl")
include("./src/run.jl")
include("./src/struct.jl")
include("./src/timestep.jl")

initialize()
run()
finish()
