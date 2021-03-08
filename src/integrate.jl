function integrate(dt)
  if time_integrator == "leap frog"
    kick_drift(0.5 * dt, dt)
    evaluate()
    kick(0.5 * dt)
  elseif time_integrator == "RK2"
    kick_drift(0.5 * dt, dt)
    evaluate()
    kick(0.5 * dt)
    evaluate()
  elseif time_integrator == "VL2"
    VL2_s0(dt)
    evaluate()
    VL2_s1(dt)
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

# van Leer predictor-corrector midpoint integrator
function VL2_s0(dt)
  Threads.@threads for i in active_particle
    x_old = sph[i].x
    p_old = sph[i].p
    U_old = sph[i].U

    # predict midpoint
    sph[i].x = x_old + p_old / sph[i].m * 0.5 * dt
    sph[i].p = p_old + sph[i].F1 * 0.5 * dt
    sph[i].U = U_old + sph[i].dU * 0.5 * dt

    # store variables
    sph[i].p_old = p_old
    sph[i].U_old = U_old

    if time_dependent_viscosity
      alpla_old = sph[i].alpha
      sph[i].alpha = alpha_old + sph[i].dalpla * 0.5 * dt
      sph[i].alpha_old = alpha_old
    end
  end
end

function VL2_s1(dt)
  Threads.@threads for i in active_particle
    p_old = sph[i].p_old
    U_old = sph[i].U_old
    x_old = sph[i].x - p_old / sph[i].m * 0.5 * dt
    sph[i].x = x_old + sph[i].p / sph[i].m * dt
    sph[i].p = p_old + sph[i].F1 * dt
    sph[i].U = U_old + sph[i].dU * dt

    if time_dependent_viscosity
      sph[i].alpha = sph[i].alpha_old + sph[i].dalpha * dt
    end
  end
end