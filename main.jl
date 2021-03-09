# simulation settings

# kernel function: "cubic spline", "Wendland C2"
kernel = "Wendland C2"

# kernel interpolation: "standard"
gradient = "standard"

# sph formulation: "vanilla ice", "Lagrangian"
formulation = "Lagrangian"

# volume element: "mass", "U"
volume_element = "mass"

# time integrator: "RK2", "VL2", "RK3", "leapfrog", "modified leapfrog", "symplectic Euler"
time_integrator = "RK3"

# output file name
outputfile = "out/out"

# The Courant-Friedrichs-Levi (CFL) number
cfl = 0.3

# time-dependent viscosity. if set to false, viscosity alpha is always 1.
time_dependent_viscosity = true
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
x_min = -1.0
x_max = 1.0
x_center = 0

# problem setup: "standard", "strong", "custom". parameters are set below.
setup = "standard"

# plot density, pressure and velocity: true, false
plot_figure = true
x_min_plot = -0.5
x_max_plot = 0.5

# number of sample points for exact Riemann solver 
Nsample_riemann = 1e3

# tolerance criterion in iterative process in exact Riemann solver and density estimate
TOL = 1e-6

# check conservation of momentum and energy.
# If this option set to true, wall particles are removed and simulated particles expand to vacuum.
debug = true

# standard shock tube (Saitoh & Makino, 2013; Springel et al., 2020)
if setup == "standard"
rho_left = 1
v_left = 0
P_left = 1
rho_right = 0.25
v_right = 0
P_right = 0.1795
# strong shock tube (Saitoh & Makino, 2013)
elseif setup == "strong"
rho_left = 1
v_left = 0
P_left = 1000
rho_right = 1
v_right = 0
P_right = 0.01
elseif setup == "custom"
# density and pressure in the left side
rho_left = 1
v_left = 0
P_left = 1000
# density and pressure in the right side
rho_right = 1
v_right = 0
P_right = 0.01
else 
  @assert 0
end

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
