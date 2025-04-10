function csbUpdate(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int; timeStep::Float64 = 0.1, 
    discreteState::Bool = false)
    n_nodes = size(update_matrix,1)
    stateVec = [-1, 1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    update_matrix = float(update_matrix)
    degreeVec = sum(abs.(update_matrix), dims = 1)
    degreeVec = 0.25.*(degreeVec.^0.5)
    for i in 1:nInit
        if discreteState
            state = float(rand(stateVec, n_nodes))
            init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        else
            state = (rand(n_nodes) - 0.5)/0.5
            init = join(["'", join(state, ","), "'"])
        end
        push!(initVec, init)
        update_matrix2 = sparse(update_matrix')
        flag = 0
        time = 0
        for j in 1:nIter
            time = time + 1
            activation = update_matrix2*state
            s1 = state .+ timeStep.*(activation./(abs.(activation) .+ 1))
            if (!(any(activation.*s1 .< 0) && any(abs.(activation) .>= degreeVec)))
                flag = 1
                break
            end
            u = rand(1:n_nodes, 1)
            state[u] = s1[u]
        end
        state = convert.(Int, sign.(state))
        fin = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        fr = frustration(state, findnz(sparse(update_matrix)))
        push!(finVec, fin)
        push!(flagVec, flag)
        push!(frustVec, fr)
        push!(timeVec, time)
    end
    states_df = DataFrame(init=initVec, 
            fin = finVec, flag = flagVec)
    frust_df = DataFrame(fin = finVec, 
        frust = frustVec)
    frust_df = unique(frust_df, :fin)
    timeData = DataFrame(fin = finVec, time = timeVec)
    timeData = groupby(timeData, :fin)
    timeData = combine(timeData, :time => avg, renamecols = false)
    frust_df = innerjoin(frust_df, timeData, on = :fin)
    # print(frust_df)
    return states_df, frust_df
end