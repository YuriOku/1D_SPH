function initialize()
  println("begin.")
  global t = t_start
  global t_output = output_interval
  global output_count = 1

  global lbox = x_max - x_min
  global sph = Array{Ptype,1}(undef, Npart)

  make_shock_tube()

  initialize_exact_riemann_solver()

  if debug
    global p_initial = 0
    global E_initial = 0
    for i = 1:length(sph)
      p_initial += sph[i].p
      E_initial += sph[i].U + 0.5 * sph[i].p^2 / sph[i].m
    end
  end

  Nthreads = Threads.nthreads()
  return println("Number of threads: $Nthreads")
end

function make_shock_tube()
  mass_left = rho_left * (center - x_min)
  mass_right = rho_right * (x_max - center)
  mass_total = mass_left + mass_right

  particle_mass = mass_total / Npart

  Nleft = ceil(mass_left / particle_mass)

  interval_left = particle_mass / rho_left
  interval_right = particle_mass / rho_right

  for i = 1:Npart
    if i < Nleft
      x = x_min + interval_left * (i - 0.5)
      rho = rho_left
      v = v_left
      P = P_left
    else
      x = center + interval_right * ((i - Nleft) - 0.5)
      rho = rho_right
      v = v_right
      P = P_right
    end
    m = particle_mass
    sph[i] = Ptype(x, v, m, rho, P)
  end

  if !debug
    global active_particle = Nngb:(Npart-Nngb)
  else
    global active_particle = 1:Npart
  end
end

function initialize_exact_riemann_solver()
  global P_s = P_star()
  global v_s = v_star(P_s)

  step_riemann = lbox / Nsample_riemann
  x_riemann = x_min:step_riemann:x_max
  rho_riemann_array = map(x -> rho_riemann(x, t_end, v_s, P_s, center), x_riemann)
  v_riemann_array = map(x -> v_riemann(x, t_end, v_s, P_s, center), x_riemann)
  P_riemann_array = map(x -> P_riemann(x, t_end, v_s, P_s, center), x_riemann)

  # y-axis limits
  global rho_max = maximum(rho_riemann_array)
  global rho_min = minimum(rho_riemann_array)
  global v_max = maximum(v_riemann_array)
  global v_min = minimum(v_riemann_array)
  global P_max = maximum(P_riemann_array)
  global P_min = minimum(P_riemann_array)
  rho_max = rho_max + 0.1 * (rho_max - rho_min)
  rho_min = rho_min - 0.1 * (rho_max - rho_min)
  v_max = v_max + 0.1 * (v_max - v_min)
  v_min = v_min - 0.1 * (v_max - v_min)
  P_max = P_max + 0.1 * (P_max - P_min)
  P_min = P_min - 0.1 * (P_max - P_min)
end