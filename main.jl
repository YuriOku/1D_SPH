# simulation settings

# kernel function: "cubic spline"
kernel = "cubic spline"

# kernel interpolation: "standard"
gradient = "standard"

# sph formulation: "vanilla ice"
formulation = "vanilla ice"

# volume element: "mass", "U"
volume_element = "U"

# output file name
outputfile = "out/out"

# The Courant-Friedrichs-Levi (CFL) number
cfl = 0.3

# upper limit of time step
dt_max = 1

# number of SPH neighbour particles
Nngb = 32

# simulation time and output timing
t_start = 0
t_end = 0.25
output_interval = 0.01

# specific heat ratio
gamma = 1.4

# number of particles
Npart = 500

# simulation box boundaries
x1_min = -1
x1_max = 1
center = 0

# density and pressure in the left side
rho_left = 1
P_left = 1

# density and pressure in the right side
rho_right = 0.25
P_right = 0.1795

# plot density, pressure and velocity: true, false
plot_figure = false

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
include("./src/kernel.jl")
include("./src/timestep.jl")
include("./src/output.jl")

initialize()
run()
finish()
