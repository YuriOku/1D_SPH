function integrate(dt)
  if time_integrator == "Leapfrog"
    kick_drift(0.5 * dt, dt)
    evaluate()
    kick(0.5 * dt)
  elseif time_integrator == "ModifiedLeapfrog"
    kick_drift(0.5 * dt, dt)
    evaluate()
    kick(0.5 * dt)
    evaluate()
  elseif time_integrator == "RK2"
    RK(0, 1, 1, dt, 1)
    evaluate()
    RK(0.5, 0.5, 0.5, dt, 2)
    evaluate()
  elseif time_integrator == "VL2"
    RK(0, 1, 0.5, dt, 1)
    evaluate()
    RK(0, 1, 1, dt, 2)
    evaluate()
  elseif time_integrator == "RK3"
    RK(0, 1, 1, dt, 1)
    evaluate()
    RK(0.25, 0.75, 0.25, dt, 2)
    evaluate()
    RK(2/3, 1/3, 2/3, dt, 3)
    evaluate()
  elseif time_integrator == "symplectic Euler"
    symplectic_euler(dt)
    evaluate()
  else
    @assert 0
  end
end


function kick_drift(dt_kick, dt_drift)
  Threads.@threads for i in active_particle
    p_old = sph[i].p
    U_old = sph[i].U

    # v(t + 0.5dt) = v(t) + a(t)*dt/2
    p_new = p_old + sph[i].F1 * dt_kick
    U_new = U_old + sph[i].dU * dt_kick

    # r(t + dt) = r(t) + v(t + 0.5dt)*dt
    sph[i].x = sph[i].x + p_new / sph[i].m * dt_drift

    # store variables for correction in the next kick
    sph[i].p_old = p_new
    sph[i].U_old = U_new

    # re-integrate for evaluation in the next step (predict)     
    sph[i].p = p_old + sph[i].F1 * dt_drift
    sph[i].U = U_old + sph[i].dU * dt_drift

    if time_dependent_viscosity
      alpha_old = sph[i].sph[i].alpha
      sph[i].alpha_old = alpha_old + sph[i].dalpha * dt_kick
      sph[i].alpha = alpha_old + sph[i].dalpha * dt_drift
    end
  end
end

function kick(dt_kick)
  Threads.@threads for i in active_particle
    # integrate again by using a(t + dt) (correct)    
    # v(t + dt) = v(t + 0.5dt) + a(t + dt)*dt/2
    sph[i].p = sph[i].p_old + sph[i].F1 * dt_kick
    sph[i].U = sph[i].U_old + sph[i].dU * dt_kick

    if time_dependent_viscosity
      sph[i].alpha = sph[i].alpha_old + sph[i].dalpha * dt_kick
    end
  end
end

# the same as kick_drift(dt, dt)
function symplectic_euler(dt)
  Threads.@threads for i in active_particle
    sph[i].p = sph[i].p + sph[i].F1 * dt
    sph[i].U = sph[i].U + sph[i].dU * dt
    if time_dependent_viscosity
      sph[i].alpha = sph[i].dalpha * dt
    end
    sph[i].x = sph[i].x + sph[i].p / sph[i].m * dt
  end
end

# Algorithm 3 of Ketcheson (2010). See also Stone et al. (2020).
function RK(g1, g2, beta, dt, stage)
  Threads.@threads for i in active_particle
    if stage == 1
      x_old = sph[i].x
      p_old = sph[i].p
      U_old = sph[i].U
      alpha_old = sph[i].alpha
      sph[i].x_old = x_old
      sph[i].p_old = p_old
      sph[i].U_old = U_old
      sph[i].alpha_old = alpha_old
    else
      x_old = sph[i].x_old
      p_old = sph[i].p_old
      U_old = sph[i].U_old
      alpha_old = sph[i].alpha_old
    end

    sph[i].x = g1 * sph[i].x + g2 * x_old + sph[i].p / sph[i].m * beta * dt
    sph[i].p = g1 * sph[i].p + g2 * p_old + sph[i].F1 * beta * dt
    sph[i].U = g1 * sph[i].U + g2 * U_old + sph[i].dU * beta * dt
    sph[i].alpha = g1 * sph[i].alpha + g2 * alpha_old + sph[i].dalpha * beta * dt
  end
end