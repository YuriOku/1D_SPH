function run()
    println("running...")

  # initial step
    evaluate()

  # main loop
    while t < t_end
        dt = timestep()

        Printf.@printf("t = %5f, dt = %.5e\n", t, dt)

        kick_drift(dt)
    
        evaluate()

        kick(dt)

        output(t, t_output)

        update_time(t, dt)
    end
end

function evaluate()
    reorder()

    Threads.@threads for i in 1:Npart
        kernel_summation(i)
    end

    Threads.@threads for i in 1:Npart
    # neighbour loop
        j = index_list[i]
        j_min = max(1, j - Nngb)
        j_max = min(Npart, j + Nngb)
        for j in j_min:j_max
            k = index_order[j]
            sph[i].F1 += force(i, k)
            sph[i].dU += dU(i, k)
            if time_dependent_viscosity
                sph[i].divv += divv(i, k)
            end
        end
        if time_dependent_viscosity
            sph[i].dalpha = dalpha_visc(i)
            sph[i].divv = 0
        end
    end
end

function kick_drift(dt)
    Threads.@threads for i in active_particle

        # v(t + 0.5dt) = v(t) + a(t)*dt/2
        sph[i].p = sph[i].p + sph[i].F1 * 0.5 * dt
        sph[i].U  = sph[i].U  + sph[i].dU * 0.5 * dt

        # r(t + dt) = r(t) + v(t + 0.5dt)*dt
        sph[i].x = sph[i].x + sph[i].p / sph[i].m * dt

        # refresh
        sph[i].F1 = 0
        sph[i].dU = 0

        if time_dependent_viscosity
            sph[i].alpha = sph[i].alpha + sph[i].dalpha * dt
            sph[i].dalpha = 0
        end
    end
end

function kick(dt)
    Threads.@threads for i in active_particle
        # v(t + dt) = v(t + 0.5dt) + a(t + dt)*dt/2
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