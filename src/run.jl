function run()
  println("running...")

  # initial step
  evaluate()

  # main loop
  while t < t_end
    dt = timestep()

    std_output(dt)

    integrate(dt)

    output(t, t_output)

    update_time(t, dt)
  end
end

function std_output(dt)
  diagnose()
  Printf.@printf("t = %5f, dt = %.5e, p_err = %.8e, E_err = %.8e, L1_err = %.8e\n", t, dt, p_error,
                   E_error, L1_error)
  log_history(dt)
end

function diagnose()
  global L1_error = measure_L1error()
  p_sum = 0
  p_sum_abs = 0
  E_sum = 0
  E_sum_abs = 0
  for i in 1:length(sph)
    p_sum += sph[i].p
    p_sum_abs += abs(sph[i].p)
    E_sum += sph[i].U + 0.5 * sph[i].p^2 / sph[i].m
    E_sum_abs += abs(sph[i].U + 0.5 * sph[i].p^2 / sph[i].m)
  end
  global p_error = abs(p_sum - p_initial) / p_sum_abs
  global E_error = abs(E_sum - E_initial) / E_sum_abs
end

function measure_L1error()
  # exact solution of Riemann problem
  step_riemann = lbox / Nsample_riemann
  x_riemann = x_min + step_riemann/2:step_riemann:x_max - step_riemann/2
  rho_riemann_array = map(x -> rho_riemann(x, t, v_s, P_s, x_center), x_riemann)

  # calculate mean density in each bin
  numbin = ceil(Int, (x_max - x_min)/L1_binsize)
  rhobar_exact = zeros(numbin)
  count_exact = zeros(numbin)
  rhobar_sim = zeros(numbin)
  count_sim = zeros(numbin)
  for i in 1:length(x_riemann)
    index = ceil(Int, (x_riemann[i] - x_min)/L1_binsize)
    rhobar_exact[index] += rho_riemann_array[i]
    count_exact[index] += 1
  end
  for i in 1:length(sph)
    index = ceil(Int, (sph[i].x - x_min)/L1_binsize)
    rhobar_sim[index] += sph[i].rho
    count_sim[index] += 1
  end
  for i in 1:numbin
    rhobar_exact[i] = count_exact[i] > 0 ? rhobar_exact[i]/count_exact[i] : 0
    rhobar_sim[i] = count_sim[i] > 0 ? rhobar_sim[i]/count_sim[i] : 0
  end

  # measure L1 error of [x_min + (1/4)lbox, x_min + (3/4)lbox]
  L1error = 0
  bin25 = ceil(Int, 0.25*lbox/L1_binsize)
  bin75 = ceil(Int, 0.75*lbox/L1_binsize)
  for i in bin25:bin75
    L1error += abs(rhobar_exact[i] - rhobar_sim[i])/(bin75 - bin25)
  end
  return L1error
end

function log_history(dt)
  open(hstfilename, "a") do io
    println(io, "$t $dt $p_error $E_error $L1_error")
  end
end