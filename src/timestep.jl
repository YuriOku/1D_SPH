function timestep()
  dt = dt_max
  for i in active_particle
    @assert sph[i].P > 0

    s = index_list[i]
    vmax = 0
    for t = (s - Nngb) : (s + Nngb)
      t = order_periodic(t)
      k = index_order[t]
      vsig = cs(i) + cs(k) - 3 * eij(i, k) * vij(i, k)
      if vsig > vmax
        vmax = vsig
      end
    end
    dt_i = cfl * sph[i].hsml / vmax
    @assert dt_i > 0
    dt = min(dt, dt_i)
  end
  return dt
end

function update_time(t, dt)
  global t += dt
end