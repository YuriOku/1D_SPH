function output(t, t_output)
  if t > t_output
    x1_array = []
    v1_array = []
    rho_array = []
    P_array = []
    output_number = string(output_count)
    outputfilename = outputfile*output_number*".txt"
    open(outputfilename, "w") do io
      for i in 1:Npart
        x1 = sph[i].x1
        v1 = sph[i].p1/sph[i].m
        rho = sph[i].rho
        P = sph[i].P
        push!(x1_array, x1)
        push!(v1_array, v1)
        push!(rho_array, rho)
        push!(P_array, P)
        println(io, "$x1 $v1 $rho $P")
      end
    end
    if plot_figure
      Plots.scatter(x1_array, rho_array)
      Plots.savefig(outputfile*"_rho"*output_number*".png")
      Plots.scatter(x1_array, P_array)
      Plots.savefig(outputfile*"_P"*output_number*".png")
      Plots.scatter(x1_array, v1_array)
      Plots.savefig(outputfile*"_v"*output_number*".png")
    end
  global t_output += output_interval
  global output_count += 1
end
end