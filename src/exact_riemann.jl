# This file contains functions to calculate an analytical solution of Riemann problem.
# 
# Reference: Toro, Eleuterio F. Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media, 2013.

function P_star()
  # initial guess
  P0 = 0.5 * (P_left + P_right)

  # Newton-Raphson iterative procedure
  CHA = 1
  while CHA > TOL
    P1 = P0 - f(P0) / dfdP(P0)
    CHA = abs(P1 - P0) / 0.5 / (P1 + P0)
    P0 = P1
    @show P0
  end
  return P0
end

function v_star(P)
  return 0.5 * (v_left + v_right) + 0.5(f_right(P) - f_left(P))
end

function region(x, t, v, P, center = 0)
  # 1: left 
  # 2: right
  # 3: left star (shock)
  # 4: right star (shock)
  # 5: left star (rarefaction)
  # 6: right star (rarefaction)
  # 7: left rarefaction fan
  # 8: right rarefaction fan

  S = (x - center) / t
  if S < v # left from contact discontinuity
    if P > P_left # shock
      if S < S_shock_left(P) # outside shock
        return 1 #"left"
      else # inside shock
        return 3 #"left star"
      end
    else # rarefaction
      if S < S_head_left() # outside
        return 1 #"left"
      elseif S < S_tail_left(P) # inside rarefaction fan
        return 7 #"left rarefaction fan"
      else
        return 5 #"left star"
      end
    end
  else # right from contact discontinuity
    if P > P_right
      if S > S_shock_right(P)
        return 2 #"right"
      else
        return 4 #"right star"
      end
    else
      if S > S_head_right()
        return 2 #"right"
      elseif S > S_tail_right(P)
        return 8 #"right rarefaction fan"
      else
        return 6 #"right star"
      end
    end
  end
end

function rho_riemann(x, t, v, P, center = 0)
  reg = region(x, t, v, P, center)
  if reg == 1 #"left"
    return rho_left
  elseif reg == 2 #"right"
    return rho_right
  elseif reg == 3 # left star (shock)
    return rho_left * (P / P_left + (gamma - 1) / (gamma + 1)) /
           ((gamma - 1) / (gamma + 1) * P / P_left + 1)
  elseif reg == 4 # right star (shock)
    return rho_right * (P / P_right + (gamma - 1) / (gamma + 1)) /
           ((gamma - 1) / (gamma + 1) * P / P_right + 1)
  elseif reg == 5 #"left star (rarefaction)"
    return rho_left * (P / P_left)^(1 / gamma)
  elseif reg == 6 #"right star (rarefaction"
    return rho_right * (P / P_right)^(1 / gamma)
  elseif reg == 7 #"left rarefaction fan"
    return rho_left *
           (
      2 / (gamma + 1) + (gamma - 1) / (gamma + 1) / aL() * (v_left - x / t)
    )^(2 / (gamma - 1))
  elseif reg == 8 #"right rarefaction fan"
    return rho_right *
           (
      2 / (gamma + 1) - (gamma - 1) / (gamma + 1) / aR() * (v_right - x / t)
    )^(2 / (gamma - 1))
  else
    @assert 0
  end
end

function v_riemann(x, t, v, P, center = 0)
  reg = region(x, t, v, P, center)
  if reg == 1 #"left"
    return v_left
  elseif reg == 2 #"right"
    return v_right
  elseif reg == 3 || reg == 4 || reg == 5 || reg == 6 #"left star" || "right star"
    return v
  elseif reg == 7 #"left rarefaction fan"
    return 2 / (gamma + 1) * (aL() + (gamma - 1) / 2 * v_left + x / t)
  elseif reg == 8 #"right rarefaction fan"
    return 2 / (gamma + 1) * (-aR() + (gamma - 1) / 2 * v_right + x / t)
  else
    @assert 0
  end
end

function P_riemann(x, t, v, P, center = 0)
  reg = region(x, t, v, P, center)
  if reg == 1 #"left"
    return P_left
  elseif reg == 2 #"right"
    return P_right
  elseif reg == 3 || reg == 4 || reg == 5 || reg == 6 #"left star" || "right star"
    return P
  elseif reg == 7 #"left rarefaction fan"
    return P_left *
           (
      2 / (gamma + 1) + (gamma - 1) / (gamma + 1) / aL() * (v_left - x / t)
    )^(2 * gamma / (gamma - 1))
  elseif reg == 8 #"right rarefaction fan"
    return P_right *
           (
      2 / (gamma + 1) - (gamma - 1) / (gamma + 1) / aR() * (v_right - x / t)
    )^(2 * gamma / (gamma - 1))
  else
    @assert 0
  end
end

function f_left(P)
  if P > P_left # shock
    f = (P - P_left) * sqrt(AL() / (P + BL()))
  else # rarefaction
    f = 2 * aL() / (gamma - 1) * ((P / P_left)^((gamma - 1) / 2 / gamma) - 1)
  end
  return f
end

function f_right(P)
  if P > P_right # shock
    f = (P - P_right) * sqrt(AR() / (P + BR()))
  else # rarefaction
    f = 2 * aR() / (gamma - 1) * ((P / P_right)^((gamma - 1) / 2 / gamma) - 1)
  end
  return f
end

function dfdP_left(P)
  if P > P_left
    dfdP = sqrt(AL() / (BL() + P)) * (1 - (P - P_left) / 2 / (BL() + P))
  else
    dfdP = 1 / rho_left / aL() * (P / P_left)^(-(gamma + 1) / 2 / gamma)
  end
  return dfdP
end

function dfdP_right(P)
  if P > P_right
    dfdP = sqrt(AR() / (BR() + P)) * (1 - (P - P_right) / 2 / (BR() + P))
  else
    dfdP = 1 / rho_right / aR() * (P / P_right)^(-(gamma + 1) / 2 / gamma)
  end
  return dfdP
end

function AL()
  return 2 / (gamma + 1) / rho_left
end
function AR()
  return 2 / (gamma + 1) / rho_right
end
function BL()
  return (gamma - 1) / (gamma + 1) * P_left
end
function BR()
  return (gamma - 1) / (gamma + 1) * P_right
end
function aL()
  return sqrt(gamma * P_left / rho_left)
end
function aR()
  return sqrt(gamma * P_right / rho_right)
end

function f(P)
  return f_left(P) + f_right(P) + v_right - v_left
end

function dfdP(P)
  return dfdP_left(P) + dfdP_right(P)
end

function S_shock_left(P)
  return v_left -
         aL() * ((gamma + 1) / 2 / gamma * P / P_left + (gamma - 1) / 2 / gamma)^(1 / 2)
end

function S_shock_right(P)
  return v_right +
         aR() * ((gamma + 1) / 2 / gamma * P / P_right + (gamma - 1) / 2 / gamma)^(1 / 2)
end

function S_head_left()
  return v_left - aL()
end
function S_tail_left(P)
  rho_s = rho_left * (P / P_left)^(1 / gamma)
  v_s = v_star(P)
  aL_s = sqrt(gamma * P / rho_s)
  return v_s - aL_s
end
function S_head_right()
  return v_right + aR()
end
function S_tail_right(P)
  rho_s = rho_right * (P / P_right)^(1 / gamma)
  v_s = v_star(P)
  aR_s = sqrt(gamma * P / rho_s)
  return v_s + aR_s
end
