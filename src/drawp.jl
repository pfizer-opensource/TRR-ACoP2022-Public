function drawp(p0::Vector{Float64}, dparam::Distribution,
    fiti::Vector{Int64}, log_flag::Bool)::Vector{Float64}
    p = Vector{Float64}(undef,size(p0)[1])
    p .= copy(p0)

    r = rand(dparam)

    for (idx, val) in enumerate(fiti)
        if log_flag
            p[val] = exp(r[idx])
        else
            p[val] = r[idx]
        end
    end
    
    return p
end

