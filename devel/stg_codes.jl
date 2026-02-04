function stg_to_df(stg, nodes)
    source = collect(keys(stg))
    target = collect(values(stg))
    source_str = [join(state, "_") for state in source]
    target_str = [join(state, "_") for state in target]
    return DataFrame([source_str, target_str], [:Source, :Target])
end

function hamming_distance(x, y)
    if (length(x) == length(y))
        return sum(x .!= y)/length(x)
    else
        return NaN
    end
end

function hamming_distance(x)
    n = length(x)
    hammings = Float64[]
    for i in 1:(n-1)
        for j in (i+1):n
            push!(hammings, hamming_distance(x[i],x[j]))
        end
    end
    return hammings
end

function compare_stgs(stg_list, nms)
    stg_combined = Dict{keytype(stg_list[1]), Vector{valtype(stg_list[1])}}()
    keys_all = unique(vcat([collect(keys(stg)) for stg in stg_list]...))
    vals_all = Vector{Vector{valtype(stg_list[1])}}()
    hammings = Vector{Vector{Float64}}()
    for key in keys_all
        vals = [get(stg, key, valtype(stg)()) for stg in stg_list]
        push!(vals_all, vals)
        push!(hammings, hamming_distance(vals))
    end
    n = length(nms)
    ids = Symbol[]
    for i in 1:(n-1)
        for j in (i+1):n 
            push!(ids, Symbol(nms[i]*"_"*nms[j]))
        end
    end
    if length(ids) == 1
        hammings = DataFrame(Hamming = reduce(vcat, hammings))
    else
        hammings = DataFrame(reduce(hcat,hammings)', ids)
    end
    vals_all = DataFrame(reduce(hcat, vals_all)', Symbol.(nms))
    keys_all = DataFrame(Source = keys_all)
    return hcat(keys_all, vals_all, hammings)
end

function compare_stgs(rules_file)
    stg_logical = Boolean.bool_to_stg(rules_file)
    topo_file = replace(rules_file, ".txt" => ".topo")
    if !isfile(topo_file)
        Boolean.boolean_to_topo(rules_file)
    end
    stg_ising = Boolean.topo_to_stg(topo_file)
    stg_nising = Boolean.topo_to_stg(topo_file; mode = :nising)
    df = compare_stgs([stg_logical, stg_ising, stg_nising], ["logical", "ising", "nising"])
    CSV.write(replace(rules_file, ".txt" => "_stg_compare.csv"), df)
end