function stateConvertOddPos(state, sOld)
    if (state < 0)
        return -1.0
    elseif (state > 0.5)
        return 1.0
    elseif (state > 0)
        return 0.5
    else
        return sOld
    end
end

function stateConvertOddNeg(state, sOld)
    if (state < -0.5)
        return -1.0
    elseif (state > 0)
        return 1.0
    elseif (state < 0)
        return -0.5
    else
        return sOld
    end
end

function oddLevels(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, negativeOdd::Bool)
    n_nodes = size(update_matrix,1)
    stateVec = [-1.0, 0.5, 1.0]
    if (negativeOdd)
        stateVec = [-1.0, -0.5, 1.0]
    end
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix = float(update_matrix)
    updOriginal = copy(update_matrix)
    for i in 1:n_nodes
        n = sum(abs.(update_matrix[:,i]))
        if n == 0
            n = 1
        end
        update_matrix[:, i] = update_matrix[:, i]/n
    end
    update_matrix2 = sparse(float(update_matrix)')
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(Int.(replace(x -> x == -1 ? 0 : abs(2*x), state)), "_"), "'"])
        flag = 0
        time = 0
        for j in 1:nIter
            time = time + 1
            s1 = update_matrix2*state
            if (negativeOdd)
                s1 = [stateConvertOddNeg(s1[i], state[i]) for i in 1:n_nodes]
            else
                s1 = [stateConvertOddPos(s1[i], state[i]) for i in 1:n_nodes]
            end
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = frustration(state, findnz(sparse(updOriginal)))
        fin = join(["'", join(Int.(replace(x -> x == -1 ? 0 : abs(2*x), state)), "_"), "'"])
        push!(frustVec, fr)
        push!(initVec, init)
        push!(finVec, fin)
        push!(flagVec, flag)
        push!(timeVec, time)       
        # push!(states_df, (init, fin, flag))
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
    