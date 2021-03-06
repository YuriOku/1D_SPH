# simulation settings

# kernel function: "cubic spline"
kernel = "cubic spline"

# kernel interpolation: "standard"
gradient = "standard"

# sph formulation: "vanilla ice"
formulation = "vanilla ice"

# volume element: "mass", "U"
volume_element = "mass"

# output file name
outputfile = "out/out"

# The Courant-Friedrichs-Levi (CFL) number
cfl = 0.3

# time-dependent viscosity
time_dependent_viscosity = true
alpha_max = 1.5
alpha_min = 0.1

# upper limit of time step
dt_max = 1

# number of SPH neighbour particles
Nngb = 8

# simulation time and output timing
t_start = 0
t_end = 0.2
output_interval = 0.01

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
rho_right = 0.3
v_right = 0
P_right = 0.3

# plot density, pressure and velocity: true, false
plot_figure = true

# y-axis limits
rho_max = 1.1
rho_min = 0.2
v_max = 1
v_min = -0.1
P_max = 1.1
P_min = 0.1

# number of sample points for exact Riemann solver 
Nsample_riemann = 1e3
# tolerance criterion in iterative process in exact Riemann solver
TOL = 1e-6

# calculation

if plot_figure
  using GR
  using Plots
  ENV["GKSwstype"] = "100"
  gr()
end
using Printf

include("./src/begin.jl")
include("./src/run.jl")
include("./src/finish.jl")
include("./src/struct.jl")
include("./src/evaluate.jl")
include("./src/force.jl")
include("./src/kernel.jl")
include("./src/timestep.jl")
include("./src/output.jl")
include("./src/exact_riemann.jl")

initialize()
run()
finish()
