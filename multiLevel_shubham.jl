function numChar(xnC, s0, levels, states, vaibhav)
    if xnC == 0
        if vaibhav
            return(0)
        else
            return(s0)
        end
    else
        sg = sign(xnC)
        xnC = abs(xnC)
        how_big_is_x = sum(xnC .> levels)
        return(sg*states[how_big_is_x])
    end
end


function stateChar(state::AbstractArray, s0, levels, states, vaibhav)
    for i in 1:length(state)
       ySc = numChar(state[i], s0[i], levels, states, vaibhav[i])
       state[i] = ySc[1]
    end
    return state
end

function shubhamFrust(state::Array{Float64,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        s = state[x]*state[y]*v
        if s < 0
            frustration = frustration + abs(s)
        end
    end
    frustration = frustration/nEdges
    return frustration
end

function stateConvert(state; nLevels = 2)
    state = Int.(nLevels*state)
    stateN = replace(x -> x < 0 ? x+nLevels : x+nLevels-1, state)
    stateN = [state[i] == 0 ? 2*nLevels : stateN[i] for i in 1:length(state)]
    state = join(["'", join(stateN, "_"), "'"])
    return state
end

function shubhamBoolean(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, discrete::Bool; nLevels::Int = 2, 
    vaibhav::Bool = false, turnOffNodes::Array{Int,1} = Int[])
    n_nodes = size(update_matrix,1)
    ls = collect(1:nLevels)
    stateVec = [-1*reverse(ls)/nLevels; ls/nLevels]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    if (vaibhav)
        if (length(turnOffNodes) == 0)
            vaibhav = [true for i in 1:n_nodes]
        else
            vaibhav = [i in turnOffNodes for i in 1:n_nodes]
        end
    else
        vaibhav = [false for i in 1:n_nodes]
    end
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix = float(update_matrix)
    updOriginal = copy(update_matrix)
    for i in 1:n_nodes
        n = length([(i,j) for j=1:n_nodes if update_matrix[j,i] == 0])        
        n = n_nodes - n
        if n == 0
            n = 1
        end
        update_matrix[:, i] = update_matrix[:, i]/n
    end
    update_matrix2 = sparse(update_matrix')
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = stateConvert(state; nLevels = nLevels)
        flag = 0
        time = 0
        ls = collect(1:nLevels)
        levels = [0; ls/nLevels]
        states = ls/nLevels
        levels = levels[1:nLevels]
        for j in 1:nIter
            time = time + 1
            st = copy(state)
            if discrete
                st = sign.(st)
            end
            s1 = stateChar(update_matrix2*st, state, levels, states, vaibhav)
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = shubhamFrust(state, findnz(sparse(updOriginal)))
        fin = stateConvert(state; nLevels = nLevels)
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
    return states_df, frust_df
end
