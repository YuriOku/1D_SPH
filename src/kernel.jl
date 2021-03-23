# In this file, kernel functions and their gradients are defined

# W(r_ij, h_k)
function W(i, j, k=i)
  h = sph[k].hsml
  invh = 1 / h
  r = dist(i, j) * invh

  if kernel == "CubicSpline"
    norm = 4 / 3 * invh
    if r > 1
      return 0
    elseif r > 0.5
      return norm * 2 * (1 - r)^3
    else
      return norm * (1 - 6 * r^2 + 6 * r^3)
    end
  elseif kernel == "WendlandC2"
    norm = 5 / 4 * invh
    if r > 1
      return 0
    else
      return norm * (1 - r)^3 * (1 + 3 * r)
    end
  else
    @assert 0
  end
end

# nabla_i W(r_ij, h_k)
function gradW(i, j, k=i)
  h = sph[k].hsml
  invh = 1 / h
  r = dist(i, j) * invh

  if gradient == "Standard"
    if kernel == "CubicSpline"
      norm = 4 / 3 * invh
      if r > 1
        gradW = 0
      elseif r > 0.5
        gradW = norm * 6 * (1 - r)^2 * (-invh)
      else
        gradW = norm * (-12 * r * invh + 18 * r^2 * invh)
      end
    elseif kernel == "WendlandC2"
      norm = 5 / 4 * invh
      if r > 1
        gradW = 0
      else
        gradW = norm * (-3 * (1 - r)^2 * (1 + 3 * r) * invh + 3 * (1 - r)^3 * invh)
      end
    else
      @assert 0
    end
    gradW = eij(i, j) * gradW
  elseif gradient == "IntegralApproach"
    gradW = - eij(i, j) * dist(i, j) * W(i, j, k) / sph[k].c11
  end
  return gradW
end

# dW(r_ij, h_i)/dh_i
function dWdh(i, j)
  h = sph[i].hsml
  invh = 1 / h
  r = dist(i, j) * invh

    if kernel == "CubicSpline"
      if r > 1
        return 0
      elseif r > 0.5
        return -8 / 3 * invh * (1 - r)^2 * (1 + r)
      else
        return -4 / 3 * invh * (1 - 18 * r^2 + 24 * r^3)
      end
    elseif kernel == "WendlandC2"
      if r > 1
        return 0
      else
        return 5 / 4 * invh^2 * (1 - r)^2 * (3 * r - 1) * (5 * r + 1)
      end
    else
      @assert 0
    end
end

# averaged for symmetry against exchange of i and j
function W_sym(i, j)
  return 0.5 * (W(i, j) + W(j, i))
end

# nabla_i (W(r_ij, h_i) + W(r_ij, h_j))/2
function gradW_sym(i, j)
  return 0.5 * (gradW(i, j) + gradW(i, j, j))
end


# This correction matrix accounts for local particle distribution.
# In this one-dimentional code, this is not a matrix but a scaler
function correction_matrix(i, j)
  return Z(j)/y(j) * dist(i, j)^2 * W(i, j)
end