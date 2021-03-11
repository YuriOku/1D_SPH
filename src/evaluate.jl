function evaluate()
  reorder()

  Threads.@threads for i = 1:Npart
    kernel_summation(i)
  end

  if gradient == "IntegralApproach"
    Threads.@threads for i = 1:Npart
      sph[i].c11 = 0
      # neighbour loop
      j = index_list[i]
      j_min = max(1, j - Nngb)
      j_max = min(Npart, j + Nngb)
      for j = j_min:j_max
        k = index_order[j]
        sph[i].c11 += correction_matrix(i, k)
      end
    end
  end

  Threads.@threads for i = 1:Npart
    # refresh
    sph[i].F1 = 0
    sph[i].dU = 0
    sph[i].dalpha = 0
    divv = 0

    # neighbour loop
    j = index_list[i]
    j_min = max(1, j - Nngb)
    j_max = min(Npart, j + Nngb)
    for j = j_min:j_max
      k = index_order[j]
      sph[i].F1 += force(i, k)
      sph[i].dU += dU(i, k)
      if time_dependent_viscosity
        divv += divergence_v(i, k)
      end
    end
    if time_dependent_viscosity
      sph[i].dalpha = dalpha_visc(i, divv)
    end
  end
end

function reorder()
  x_array = []
  for i = 1:Npart
    push!(x_array, sph[i].x)
  end
  global index_order = sortperm(x_array)
  global index_list = zeros(Int, Npart)
  for i = 1:Npart
    j = index_order[i]
    index_list[j] = i
  end
end

function kernel_summation(i)
  # search neighbour, adjust smoothing length and estimate density
  evaluate_hsml(i)
  density(i, volume_element)
  # EoS
  eos(i, volume_element)
end

function evaluate_hsml(i)
  h0 = eta_hsml * (sph[i].m / sph[i].rho)
  CHA = 1
  iter = 0
  while CHA > TOL
    sph[i].hsml = h0
    density(i, "mass")
    h1 = h0 - g(i) / dgdh(i)
    CHA = abs(h0 - h1) / 2 / (h0 + h1)
    h0 = h1
    iter += 1
    if iter > 1e4
      sph[i].hsml = h0
      println("iteration loop of hsml evaluation exceeds 1e4.")
      @printf("i = %d, h0 = %.7f, CHA = %.7f, g(i) = %.7f, rho = %.7f\n", i, h0, CHA, g(i), density(i, "mass"))
      break
    end
  end
end

function g(i)
  return sph[i].hsml - eta_hsml * sph[i].m / sph[i].rho
end

function dgdh(i)
  drhodh = 0
  sph[i].dydh = 0
  j = index_list[i]
  j_min = max(1, j - Nngb)
  j_max = min(Npart, j + Nngb)
  for j = j_min:j_max
    k = index_order[j]
    drhodh += sph[k].m * dWdh(i, k)
    sph[i].dydh += Z(k) * dWdh(i, k)
  end
  dgdh = 1 + eta_hsml * sph[i].m / sph[i].rho^2 * drhodh
  return dgdh
end

function density(i::Int, Z_input)
  y = 0
  j = index_list[i]
  j_min = max(1, j - Nngb)
  j_max = min(Npart, j + Nngb)
  for j = j_min:j_max
    k = index_order[j]
    if Z_input == "mass"
      y += sph[k].m * W(i, k)
    elseif Z_input == "U"
      y += sph[k].U * (gamma - 1) * W(i, k)
    else
      @assert 0
    end
  end

  if Z_input == "mass"
    sph[i].rho = y
  elseif Z_input == "U"
    sph[i].P = y
  else
    @assert 0
  end
end

function eos(i::Int, Z_input)
  if sph[i].U < 0
    @show i, sph[i].x, sph[i].p / sph[i].m
  end
  @assert sph[i].U > 0
  if Z_input == "mass"
    sph[i].P = (gamma - 1) * sph[i].rho / sph[i].m * sph[i].U
  elseif Z_input == "U"
    sph[i].rho = sph[i].P / (gamma - 1) * sph[i].m / sph[i].U
  else
    @assert 0
  end
end

