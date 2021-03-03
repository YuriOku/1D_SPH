function initialize()
  println("begin.")
  global t = t_start
  global t_output = output_interval
  global output_count = 0

  global lbox = x1_max - x1_min
  global sph = Array{Ptype, 1}(undef, Npart)

  global p_initial = 0
  global U_initial = P_left/(gamma - 1)*(center - x1_min) + P_right/(gamma - 1)*(x1_max - center)

  make_shock_tube()

  Nthreads = Threads.nthreads()
  println("Number of threads: $Nthreads")
end

function make_shock_tube()
  mass_left = rho_left*(center - x1_min)
  mass_right = rho_right*(x1_max - center)
  mass_total = mass_left + mass_right

  particle_mass = mass_total / Npart

  Nleft = ceil(mass_left/particle_mass)

  interval_left = particle_mass / rho_left
  interval_right = particle_mass / rho_right

  for i in 1:Npart
    if i < Nleft
      x1 = x1_min + interval_left*(i - 0.5)
      rho = rho_left
      P = P_left
    else
      x1 = center + interval_right*((i - Nleft) - 0.5)
      rho = rho_right
      P = P_right
    end
    v1 = 0
    m = particle_mass
    sph[i] = Ptype(x1, v1, m, rho, P)
  end

  global active_particle = Nngb:Npart - Nngb
end