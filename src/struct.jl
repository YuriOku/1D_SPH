# In this file, data structure of particle and basic functions are defined.

mutable struct Ptype
  x::Float64
  p::Float64
  m::Float64
  rho::Float64
  U::Float64
  P::Float64

  F1::Float64
  dU::Float64
  hsml::Float64

  divv::Float64
  alpha::Float64
  dalpha::Float64

  function Ptype(x, v, m, rho, P)
    p = m*v
    U = P/(gamma - 1)*m/rho
    return new(x, p, m, rho, U, P, 0, 0, 0, 0, alpha_min, 0)
  end
end

function Z(i::Int)
  if volume_element == "mass"
    return sph[i].m
  elseif volume_element == "U"
    return sph[i].U
  end
end
function y(i::Int)
  if volume_element == "mass"
    return sph[i].rho
  elseif volume_element == "U"
    return sph[i].P/(gamma - 1)
  end
end

function v(i::Int)
  return sph[i].p/sph[i].m
end

function dist(i::Int, j::Int)
  rij = abs(sph[i].x - sph[j].x)
  return rij 
end

function eij(i::Int, j::Int)
  rij = abs(sph[i].x - sph[j].x)
  eij = (sph[i].x - sph[j].x)/rij

  return eij
end

function vij(i::Int, j::Int)
  vij = v(i) - v(j)

  return vij
end

function cs(i)
  cs = sqrt(gamma * sph[i].P / sph[i].rho)
  return cs
end