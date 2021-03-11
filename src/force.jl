function force(i::Int, j::Int)
  if i == j
    return 0
  end

  if formulation == "VanillaIce"
    force_hydro =
      -Z(i) * Z(j) / y(i) / y(j) * (sph[i].P + sph[j].P) * gradW_sym(i, j)
  elseif formulation == "Lagrangian"
    force_hydro =
      -Z(i) * Z(j) * (
        grad_h(i) * sph[i].P / y(i)^2 * gradW(i, j) +
        grad_h(j) * sph[j].P / y(j)^2 * gradW(i, j, j)
      )
  else
    @assert 0
  end
  force_visc = -sph[i].m * sph[j].m * PI_visc(i, j) * gradW_sym(i, j)
  return force_hydro + force_visc
end

function dU(i::Int, j::Int)
  if i == j
    return 0
  end

  if formulation == "VanillaIce"
    dU_hydro =
      sph[i].P * Z(i) * Z(j) / y(i) / y(j) * vij(i, j) * gradW_sym(i, j)
  elseif formulation == "Lagrangian"
    dU_hydro =
      grad_h(i) * sph[i].P * Z(i) * Z(j) / y(i)^2 * vij(i, j) * gradW(i, j)
  else
    @assert 0
  end
  dU_visc =
    0.5 * sph[i].m * sph[j].m * PI_visc(i, j) * vij(i, j) * gradW_sym(i, j)
  return dU_hydro + dU_visc
end

function grad_h(i)
  invf = 1 + sph[i].hsml / y(i) * sph[i].dydh
  return 1/invf
end

function PI_visc(i, j)
  if eij(i, j) * vij(i, j) > 0
    return 0
  else
    if !time_dependent_viscosity
      alpha = 1
      beta = 2
    else
      alpha = (sph[i].alpha + sph[j].alpha) / 2
      beta = 2 * alpha
    end
    epsilon = 0.01
    hbar = (sph[i].hsml + sph[j].hsml) / 4
    cbar = (cs(i) + cs(j)) / 2
    rhobar = (sph[i].rho + sph[j].rho) / 2
    mu = (hbar * dist(i, j) * eij(i, j) * vij(i, j)) / (dist(i, j)^2 + epsilon * hbar^2)
    PI_visc = (-alpha * cbar * mu + beta * mu^2) / rhobar
    return PI_visc
  end
end

function divergence_v(i::Int, j::Int)
  if i == j
    return 0
  else
    return -Z(j) / y(j) * vij(i, j) * eij(i, j) * gradW_sym(i, j)
  end
end

function dalpha_visc(i, divv)
  S = max(-divv * (alpha_max - sph[i].alpha), 0)
  tau = sph[i].hsml / 0.2 / cs(i)
  dalpha_visc = -(sph[i].alpha - alpha_min) / tau + S
  return dalpha_visc
end