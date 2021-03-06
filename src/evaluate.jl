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
    density(i)
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

function density(i::Int)
    dist_list = []
    
    j = index_list[i]
    j_min = max(1, j - Nngb)
    j_max = min(Npart, j + Nngb)
    for j in j_min:j_max
        k = index_order[j]
        rik = dist(i, k)
        push!(dist_list, rik)
    end
    # for j in 1:Npart
    #     rij = dist(i,j)
    #     push!(dist_list, rij)
    # end
    sph[i].hsml = sort(dist_list)[Nngb]

    if volume_element == "mass"
        rho = 0    
        j = index_list[i]
        j_min = max(1, j - Nngb)
        j_max = min(Npart, j + Nngb)
        for j in j_min:j_max
            k = index_order[j]
            rho += sph[k].m * W(i, k)
        end
        sph[i].rho = rho
    elseif volume_element == "U"
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