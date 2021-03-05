function run()
    println("running...")

  # initial step
    evaluate()

  # main loop
    while t < t_end
        dt = timestep()

        diagnose()

        Printf.@printf("t = %5f, dt = %.5f, p_err = %.5e, U_err = %.5e\n", t, dt, p_error, U_error)

        kick_drift(dt)
    
        evaluate()

        kick(dt)

        output(t, t_output)

        update_time(t, dt)
    end
end

function evaluate()
    for i in 1:Npart
        kernel_summation(i)
    end

    Threads.@threads for i in 1:Npart
    # neighbour loop
        for j in 1:Npart
            sph[i].F1 += force(i, j)
            sph[i].dU += dU(i, j)

            if force(i, j) + force(j,i) > 1e-6
                @show i, j, force(i, j)/force(j,i)
            end
        end
    end
end

function kick_drift(dt)
    for i in active_particle

    # v(t + 0.5dt) = v(t) + a(t)*dt/2
        sph[i].p = sph[i].p + sph[i].F1 * 0.5 * dt
        sph[i].U  = sph[i].U  + sph[i].dU * 0.5 * dt

        # r(t + dt) = r(t) + v(t + 0.5dt)*dt
        sph[i].x = sph[i].x + sph[i].p / sph[i].m * dt

    # refresh
        sph[i].F1 = 0
        sph[i].dU = 0
    end
end

function kick(dt)
    for i in active_particle
    # v(t + dt) = v(t + 0.5dt) + a(t + 1)*dt/2
        sph[i].p = sph[i].p + sph[i].F1 * 0.5 * dt
        sph[i].U  = sph[i].U  + sph[i].dU * 0.5 * dt
    end
end

function diagnose()
  p_sum = 0
  p_sum_abs = 0
  U_sum = 0
  U_sum_abs = 0
  for i in 1:length(sph)
    p_sum += sph[i].p
    p_sum_abs += abs(sph[i].p)
    U_sum += sph[i].U
    U_sum_abs += abs(sph[i].U)
  end
  global p_error = (p_sum - p_initial)/p_sum_abs
  global U_error = (U_sum - U_initial)/U_sum_abs
end