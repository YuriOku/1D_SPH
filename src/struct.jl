# In this file, data structure of particle and basic functions are defined.

mutable struct Ptype
  x1::Float64
  p1::Float64
  m::Float64
  rho::Float64
  U::Float64
  P::Float64

  F1::Float64
  dU::Float64
  hsml::Float64

  list::Array{Int, 1}

  function Ptype(x1, v1, m, rho, P)
    p1 = m*v1
    U = P/(gamma - 1)*m/rho
    return new(x1, p1, m, rho, U, P, 0, 0, 0, [])
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

function v1(i::Int)
  return sph[i].p1/sph[i].m
end

function dist(i::Int, j::Int)
  rij = abs(sph[i].x1 - sph[j].x1)
  return rij 
end

function eij(i::Int, j::Int)
  rij = abs(sph[i].x1 - sph[j].x1)
  eij = (sph[i].x1 - sph[j].x1)/rij

  return eij
end

function vij(i::Int, j::Int)
  vij = v1(i) - v1(j)

  return vij
end

function cs(i)
  cs = sqrt(gamma * sph[i].P / sph[i].rho)
  return cs
end