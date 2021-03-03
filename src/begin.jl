function initialize()
  println("begin.")
  global t = t_start
  global t_output = output_interval
  global output_count = 0

  linecount = 0
  open(initcond, "r") do io
    for line in eachline(io)
      linecount += 1
    end
  end

  global Npart = linecount
  global lbox = x1_max - x1_min
  global sph = Array{Ptype, 1}(undef, Npart)

  linecount = 0
  x1_array = []
  open(initcond, "r") do io
    for line in eachline(io)
      linecount += 1
      words = split(line)
      x1 = parse(Float64, words[1])
      v1 = parse(Float64, words[2])
      m = parse(Float64, words[3])
      rho = parse(Float64, words[4])
      P = parse(Float64, words[5])
      sph[linecount] = Ptype(x1, v1, m, rho, P)
      push!(x1_array, x1)
    end
  end

  if boundary == "fixed"
    global active_particle = sortperm(x1_array)[Nngb: end-Nngb]
  end

end