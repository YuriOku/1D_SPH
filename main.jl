# simulation settings
kernel = "cubic spline"
gradient = "standard"
formulation = "vanilla ice"
volume_element = "mass"
initcond = "initial.txt"
outputfile = "out/out"
boundary = "fixed"
cfl = 0.1
dt_max = 1
Nngb = 32
t_start = 0
t_end = 20
output_interval = 1e-6
gamma = 1.4
x1_min = -0.5
x1_max = 0.5
plot_figure = true

# calculation

if plot_figure
  using GR
  using Plots
  ENV["GKSwstype"] = "100"
  gr()
end

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
