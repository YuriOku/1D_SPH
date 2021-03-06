function initialize()
  println("begin.")
  global t = t_start
  global t_output = output_interval
  global output_count = 0

  global lbox = x_max - x_min
  global sph = Array{Ptype, 1}(undef, Npart)

  global p_initial = 0
  global U_initial = P_left/(gamma - 1)*(center - x_min) + P_right/(gamma - 1)*(x_max - center)

  make_shock_tube()

  initialize_exact_riemann_solver()

  Nthreads = Threads.nthreads()
  println("Number of threads: $Nthreads")
end

function make_shock_tube()
  mass_left = rho_left*(center - x_min)
  mass_right = rho_right*(x_max - center)
  mass_total = mass_left + mass_right

  particle_mass = mass_total / Npart

  Nleft = ceil(mass_left/particle_mass)

  interval_left = particle_mass / rho_left
  interval_right = particle_mass / rho_right

  for i in 1:Npart
    if i < Nleft
      x = x_min + interval_left*(i - 0.5)
      rho = rho_left
      v = v_left
      P = P_left
    else
      x = center + interval_right*((i - Nleft) - 0.5)
      rho = rho_right
      v = v_right
      P = P_right
    end
    m = particle_mass
    sph[i] = Ptype(x, v, m, rho, P)
  end

  global active_particle = Nngb:Npart - Nngb
end

function initialize_exact_riemann_solver()
  global P_s = P_star()
  global v_s = v_star(P_s)
end