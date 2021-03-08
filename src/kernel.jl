# In this file, kernel functions and their gradients are defined

function W(i, j)
  h = sph[i].hsml
  invh = 1 / h
  r = dist(i, j) * invh

  if kernel == "cubic spline"
    norm = 4 / 3 * invh
    if r > 1
      return 0
    elseif r > 0.5
      return norm * 2 * (1 - r)^3
    else
      return norm * (1 - 6 * r^2 + 6 * r^3)
    end
  elseif kernel == "Wendland C2"
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


function gradW(i, j)
  h = sph[i].hsml
  invh = 1 / h
  r = dist(i, j) * invh

  if gradient == "standard"
    if kernel == "cubic spline"
      norm = 4 / 3 * invh
      if r > 1
        return 0
      elseif r > 0.5
        return norm * 6 * (1 - r)^2 * (-invh)
      else
        return norm * (-12 * r * invh + 18 * r^2 * invh)
      end
    elseif kernel == "Wendland C2"
      norm = 5 / 4 * invh
      if r > 1
        return 0
      else
        return norm * (-3 * (1 - r)^2 * (1 + 3 * r) * invh + 3 * (1 - r)^3 * invh)
      end
    else
      @assert 0
    end
  end
end

function dWdh(i, j)
  h = sph[i].hsml
  invh = 1 / h
  r = dist(i, j) * invh

  if gradient == "standard"
    if kernel == "cubic spline"
      norm = 4 / 3 * invh
      if r > 1
        return 0
      elseif r > 0.5
        return norm * 6 * (1 - r)^2 * (r * invh)
      else
        return norm * (-12 * r * (-r * invh) + 18 * r^2 * (-r * invh))
      end
    elseif kernel == "Wendland C2"
      norm = 5 / 4 * invh
      if r > 1
        return 0
      else
        return norm *
               (-3 * (1 - r)^2 * (1 + 3 * r) * (-r * invh) + 3 * (1 - r)^3 * (-r * invh))
      end
    else
      @assert 0
    end
  end
end

function W_tilde(i, j)
  return 0.5 * (W(i, j) + W(j, i))
end

function gradW_tilde(i, j)
  return 0.5 * (gradW(i, j) + gradW(j, i))
end
