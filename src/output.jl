function output(t, t_output)
  if t > t_output
    x_array = []
    v_array = []
    rho_array = []
    P_array = []
    outputfilename = @sprintf("%s_%03d.txt", outputfile, output_count)

    open(outputfilename, "w") do io
      println(io, "### x, v, rho, P")
      for i in 1:Npart
        x = sph[i].x
        v = sph[i].p / sph[i].m
        rho = sph[i].rho
        P = sph[i].P
        hsml = grad_h(i)
        push!(x_array, x)
        push!(v_array, v)
        push!(rho_array, rho)
        push!(P_array, P)
        println(io, "$x $v $rho $P $hsml")
      end
    end

    if plot_figure
      step_riemann = lbox / Nsample_riemann
      x_riemann = x_min:step_riemann:x_max
      rho_riemann_array = map(x -> rho_riemann(x, t, v_s, P_s, x_center), x_riemann)
      v_riemann_array = map(x -> v_riemann(x, t, v_s, P_s, x_center), x_riemann)
      P_riemann_array = map(x -> P_riemann(x, t, v_s, P_s, x_center), x_riemann)
      settings = "$formulation, $gradient gradient\n$time_integrator, $kernel, η=$eta_hsml\nvolume element: $volume_element\ntime-dependent art. visc.: $time_dependent_viscosity\n"

      plt_rho = Plots.plot(; ylabel="rho", title=settings, titlefont=font(12),
                           ylims=(rho_min, rho_max), xlims=(x_min_plot, x_max_plot))
      plt_rho = Plots.plot!(x_riemann, rho_riemann_array; label="exact")
      plt_rho = Plots.scatter!(x_array, rho_array; markersize=2, label="sim")

      plt_P = Plots.plot(; ylabel="P", xlabel="x", ylims=(P_min, P_max),
                         xlims=(x_min_plot, x_max_plot))
      plt_P = Plots.plot!(x_riemann, P_riemann_array; label="exact")
      plt_P = Plots.scatter!(x_array, P_array; markersize=2, label="sim")

      plt_v = Plots.plot(; ylabel="v", ylims=(v_min, v_max), xlims=(x_min_plot, x_max_plot))
      plt_v = Plots.plot!(x_riemann, v_riemann_array; label="exact")
      plt_v = Plots.scatter!(x_array, v_array; markersize=2, label="sim")

      Plots.plot(plt_rho, plt_v, plt_P; layout=(3, 1), size=(400, 600))
      figname = @sprintf("%s_%s_%s_%s_%s_%s_%s_fig_%03d.png", outputfile, formulation, gradient, time_integrator, kernel, volume_element, time_dependent_viscosity, output_count)
      Plots.savefig(figname)
    end
    global t_output += output_interval
    global output_count += 1
  end
end