function run()
    println("running...")

  # initial step
    evaluate()

  # main loop
    while t < t_end
        dt = timestep()
        @show t, dt

        diagnose()

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

    for i in active_particle
    # neighbour loop
        for j in 1:Npart
            sph[i].F1 += force(i, j)
            sph[j].dU += dU(i, j)
        end
    end
end

function kick_drift(dt)
    for i in active_particle

    # v(t + 0.5dt) = v(t) + a(t)*dt/2
        sph[i].p1 = sph[i].p1 + sph[i].F1 * 0.5 * dt
        sph[i].U  = sph[i].U  + sph[i].dU * 0.5 * dt

    # r(t + dt) = r(t) + v(t + 0.5dt)*dt
        sph[i].x1 = sph[i].x1 + sph[i].p1 / sph[i].m * dt

    # periodic_boundary(i)

    # refresh
        sph[i].F1 = 0
        sph[i].dU = 0
    end
end

function kick(dt)
    for i in active_particle
    # v(t + dt) = v(t + 0.5dt) + a(t + 1)*dt/2
        sph[i].p1 = sph[i].p1 + sph[i].F1 * 0.5 * dt
        sph[i].U  = sph[i].U  + sph[i].dU * 0.5 * dt
    end
end

function periodic_boundary(i)
    if sph[i].x1 < x1_min
        sph[i].x1 += lbox
    elseif sph[i].x1 > x1_max
    sph[i].x1 -= lbox
    end
end

function diagnose()
  p_sum = 0
  U_sum = 0
  for i in 1:length(sph)
    p_sum += sph[i].p1
    U_sum += sph[i].U
  end
  @show p_sum, U_sum
end