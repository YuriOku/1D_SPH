function reorder()
    x_array = []
    for i in 1:Npart
        push!(x_array, sph[i].x)
    end
    global index_order = sortperm(x_array)
    global index_list = zeros(Int, Npart)
    for i in 1:Npart
        j = index_order[i]
        index_list[j] = i
    end
end

function kernel_summation(i)
    # search neighbour, adjust smoothing length and estimate density
    evaluate_hsml(i)
    density(i, volume_element)
    # EoS
    if volume_element == "mass"
        sph[i].P = eos_P(i)
    elseif volume_element == "U"
        sph[i].rho = eos_rho(i)
    end
end

function eos_P(i::Int)
    if sph[i].U < 0
        @show i, sph[i].x, sph[i].p / sph[i].m
    end
    @assert sph[i].U > 0
    P = (gamma - 1) * sph[i].rho / sph[i].m * sph[i].U
    return P
end

function eos_rho(i::Int)
    if sph[i].U < 0
        @show i, sph[i].x, sph[i].p / sph[i].m
    end
    @assert sph[i].U > 0
    rho = sph[i].P / (gamma - 1) * sph[i].m / sph[i].U
    return rho
end

function density(i::Int, volume)
    if volume == "mass"
        rho = 0    
        j = index_list[i]
        j_min = max(1, j - Nngb)
        j_max = min(Npart, j + Nngb)
        for j in j_min:j_max
            k = index_order[j]
            rho += sph[k].m * W(i, k)
        end
        sph[i].rho = rho
    elseif volume == "U"
        P = 0
        j = index_list[i]
        j_min = max(1, j - Nngb)
        j_max = min(Npart, j + Nngb)
        for j in j_min:j_max
            k = index_order[j]
            P += sph[k].U * (gamma - 1) * W(i, k)
        end
        sph[i].P = P
    end     
end

function evaluate_hsml(i)
    h0 = eta_hsml * (sph[i].m/sph[i].rho)
    CHA = 1
    iter = 0
    while CHA > TOL
        sph[i].hsml = h0
        density(i, "mass")
        h1 = h0 - g(i)/dgdh(i)
        CHA = abs(h0 - h1)/2/(h0 + h1)
        h0 = h1
        # h0 = 0.1*(h1 - h0) + h0
        iter += 1
        if iter > 1e4
            sph[i].hsml = h0
            println("iteration loop of density evaluation exceeds 1e4. h0 = $h0")
            break
        end
    end
end

function g(i)
    return sph[i].hsml - eta_hsml*sph[i].m/sph[i].rho
end

function dgdh(i)
    dgdh = 1
    j = index_list[i]
    j_min = max(1, j - Nngb)
    j_max = min(Npart, j + Nngb)
    for j in j_min:j_max
        k = index_order[j]
        dgdh += eta_hsml*sph[i].m/sph[i].rho^2 * sph[k].m * dWdh(i,k)
    end
    return dgdh
end