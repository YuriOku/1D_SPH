function timestep()
  dt = dt_max
  for i in active_particle
    @assert sph[i].P > 0
    for j in sph[i].list
      vsig = cs(i) + cs(j) - 3*eij(i,j)*vij(i,j)
      dt_ij = cfl * sph[i].hsml/vsig

      if !isnan(dt_ij) && dt_ij > 0
        dt = min(dt, dt_ij)
      end
    end
  end
  return dt
end

function update_time(t, dt)
  global t += dt
end