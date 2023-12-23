function numChar(x, s0, levels, states, vaibhav)
    if x == 0
        if vaibhav
            y = 0
        else
            y = s0
        end
    else
        compares = sum(x .>= levels)
        if compares == 0
            compares = 1
        end
        y = states[compares]
    end
end


function stateChar(state::AbstractArray, s0, levels, states, vaibhav)
    for i in 1:length(state)
        x = state[i]
        y = numChar(x, s0[i], levels, states, vaibhav)
        state[i] = y
    end
    return state
end

function shubhamFrust(state::Array{Float64,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        s = state[x]*state[y]*v
        # if s<0
        #     if state[x] == state[y]
        #         frustration = frustration + 1
        #     else
        #         frustration = frustration + 0.5
        #     end
        # end
        if s < 0
            frustration = frustration + abs(state[x]*state[y])
        end
    end
    frustration = frustration/nEdges
    return frustration
end

function stateConvert(state; nLevels = 2)
    state = Int.(nLevels*state)
    # if nLevels > 5, add "_" between the states
    # if nLevels > 5
    state = join(["'", join(replace(x -> x < 0 ? x+nLevels : x+nLevels-1, state), "_"), "'"])
    # end
    # state = join(["'", join(replace(x -> x < 0 ? x+nLevels : x+nLevels-1, state)), "'"])
    return state
end

function shubhamBoolean(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, discrete::Bool; nLevels::Int = 2, 
    vaibhav::Bool = false)
    n_nodes = size(update_matrix,1)
    ls = collect(1:nLevels)
    stateVec = [-1*reverse(ls)/nLevels; ls/nLevels]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix = float(update_matrix)
    for i in 1:n_nodes
        n = length([(i,j) for j=1:n_nodes if update_matrix[j,i] == 0])        
        n = n_nodes - n
        if n == 0
            n = 1
        end
        update_matrix[:, i] = update_matrix[:, i]/n
    end
    # update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    # nzId = enumerate(findall(update_matrix.!=0))
    # minVal = minimum([abs(update_matrix[j]) for (i,j) in nzId])
    # update_matrix2 = update_matrix + Matrix(I, n_nodes, n_nodes)*(minVal/2)
    update_matrix2 = sparse(update_matrix')
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = stateConvert(state; nLevels = nLevels)
        flag = 0
        time = 0
        ls = collect(1:nLevels)
        levels = [-1*reverse(ls)/nLevels; 0; ls/nLevels]
        states = [-1*reverse(ls)/nLevels; ls/nLevels]
        levels = levels[1:(length(levels) - 1)]
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
        fr = shubhamFrust(state, findnz(sparse(update_matrix)))
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
