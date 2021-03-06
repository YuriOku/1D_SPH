function force(i::Int, j::Int)
  if i == j
      return 0
  end
  if formulation == "vanilla ice"
      force_hydro = - Z(i) * Z(j) / y(i) / y(j) * (sph[i].P + sph[j].P) * eij(i, j) * gradW_tilde(i, j)
  end
  force_visc = - sph[i].m * sph[j].m * PI_visc(i, j) * eij(i, j) * gradW_tilde(i, j)
  return force_hydro + force_visc
end

function dU(i::Int, j::Int)
  if i == j
      return 0
  end
  if formulation == "vanilla ice"
      dU_hydro = sph[i].P * Z(i) * Z(j) / y(i) / y(j) * vij(i, j) * eij(i, j) * gradW_tilde(i, j)
  end
  dU_visc = 0.5 * sph[i].m * sph[j].m * PI_visc(i, j) * vij(i, j) * eij(i, j) * gradW_tilde(i, j)
  return dU_hydro + dU_visc
end

function PI_visc(i, j)
  if eij(i, j) * vij(i, j) > 0
      return 0
  else
      if !time_dependent_viscosity
          alpha = 1
          beta = 2
      else
          alpha = sph[i].alpha
          beta = 2 * alpha
      end
      epsilon = 0.01
      hbar = (sph[i].hsml + sph[j].hsml) / 2
      cbar = (cs(i) + cs(j)) / 2
      rhobar = (sph[i].rho + sph[j].rho) / 2
      mu = (hbar * dist(i, j) * eij(i, j) * vij(i, j)) / (dist(i, j)^2 + epsilon * hbar^2)
      PI_visc = (-alpha * cbar * mu + beta * mu^2) / rhobar
      return PI_visc
  end
end

function divv(i::Int, j::Int)
  if i == j
    return 0
  else
    return - Z(j)/y(j) * vij(i,j) * eij(i,j)*gradW_tilde(i,j) 
  end
end

function dalpha_visc(i)
  divv = sph[i].divv
  S = max(-divv*(alpha_max - sph[i].alpha), 0)
  tau = sph[i].hsml / 0.1 / cs(i)
  dalpha_visc = - (sph[i].alpha - alpha_min)/tau + S
  return dalpha_visc
end