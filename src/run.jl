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
  if !debug
    Printf.@printf("t = %5f, dt = %.5e\n", t, dt)
  else
    diagnose()
    Printf.@printf("t = %5f, dt = %.5e, p_err = %.8e, E_err = %.8e\n", t, dt, p_error,
                   E_error)
  end
end

function diagnose()
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
  return global E_error = abs(E_sum - E_initial) / E_sum_abs
end
