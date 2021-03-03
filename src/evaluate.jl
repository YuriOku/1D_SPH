function force(i::Int, j::Int)
    if i == j
        return 0
    end
    if formulation == "vanilla ice"
        force_hydro = - Z(i) * Z(j) / y(i) / y(j) * (sph[i].P + sph[j].P) * eij(i, j) * gradW_tilde(i, j)
    end
    force_visc = - sph[i].m * sph[j].m * PI_visc(i, j) * eij(i,j) * gradW_tilde(i, j)
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
        alpha = 1
        beta = 2
        epsilon = 0.01
        hbar = (sph[i].hsml + sph[j].hsml) / 2
        cbar = (cs(i) + cs(j)) / 2
        rhobar = (sph[i].rho + sph[j].rho) / 2
        mu = (hbar * dist(i, j) * eij(i, j) * vij(i, j)) / (dist(i, j)^2 + epsilon * hbar^2)
        PI_visc = (-alpha * cbar * mu + beta * mu^2) / rhobar
        return PI_visc
    end
end

function eos_P(i::Int)
    if sph[i].U < 0
        @show i, sph[i].x1, sph[i].p1 / sph[i].m
    end
    @assert sph[i].U > 0
    P = (gamma - 1) * sph[i].rho / sph[i].m * sph[i].U
    return P
end

function eos_rho(i::Int)
  if sph[i].U < 0
      @show i, sph[i].x1, sph[i].p1 / sph[i].m
  end
  @assert sph[i].U > 0
  rho = sph[i].P/(gamma - 1)*sph[i].m/sph[i].U
  return rho
end

function density(i::Int)
    dist_list = []
    for j in 1:length(sph)
        rij = dist(i, j)
        push!(dist_list, rij)
    end
    sph[i].list = sortperm(dist_list)[1:Nngb]
    hsml = dist(i, sph[i].list[end])
    sph[i].hsml = hsml

    if volume_element == "mass"
      rho = 0
      for j in sph[i].list
          rho += sph[j].m * W(i, j)
      end
      sph[i].rho = rho
    elseif volume_element == "U"
      P = 0
      for j in sph[i].list
          P += sph[j].U * (gamma - 1) * W(i, j)
      end
      sph[i].P = P
    end     
end

function kernel_summation(i)
    # search neighbour, adjust smoothing length and estimate density
    density(i)
    # EoS
    if volume_element == "mass"
      sph[i].P = eos_P(i)
    elseif volume_element == "U"
      sph[i].rho = eos_rho(i)
    end
end