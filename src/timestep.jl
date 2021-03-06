function timestep()
  dt = dt_max
  for i in active_particle
    @assert sph[i].P > 0
    j = index_list[i]
    j_min = max(1, j - Nngb)
    j_max = min(Npart, j + Nngb)
    vmax = 0
    for j in j_min:j_max
      k = index_order[j]
      vsig = cs(i) + cs(k) - 3*eij(i,k)*vij(i,k)
      if vsig > vmax
        vmax = vsig
      end
    end
    dt_i = cfl * 2 * sph[i].hsml / vmax
    @assert dt_i > 0
    dt = min(dt, dt_i)
  end
  return dt
end

function update_time(t, dt)
  global t += dt
end