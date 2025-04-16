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


function stateChar(state::Vector{Float64}, s0, levels, states, vaibhav)
    for i in 1:length(state)
       ySc = numChar(state[i], s0[i], levels[i], states[i], vaibhav[i])
       state[i] = ySc[1]
    end
    return state
end

function stateChar(state::Float64, s0, levels, states, vaibhav)
    ySc = numChar(state, s0, levels, states, vaibhav)
    return ySc[1]
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

function stateConvert(state, nLevels::Int)
    state = Int.(nLevels*state)
    stateN = replace(x -> x < 0 ? x+nLevels : x+nLevels-1, state)
    stateN = [state[i] == 0 ? 2*nLevels : stateN[i] for i in eachindex(state)]
    state = join(["'", join(stateN, "_"), "'"])
    return state
end

function stateConvert(state, nLevels::Array{Int,1})
    state = Int.(state.*nLevels)
    for i in eachindex(state)
        if state[i] < 0
            state[i] = state[i] + nLevels[i]
        elseif state[i] > 0
            state[i] = state[i] + nLevels[i] - 1
        else
            state[i] = 2*nLevels[i]
        end
    end
    state = join(state, "_")
    # end
    # state = join(["'", join(replace(x -> x < 0 ? x+nLevels : x+nLevels-1, state)), "'"])
    return state
end

function getStateVec(nLevels::Int = 2)
    ls = collect(1:nLevels)
    stateVec = float([-1*reverse(ls)/nLevels; ls/nLevels])
    return stateVec
end

function getLevels(nLevels::Int = 2)
    ls = collect(1:nLevels)
    levels = float([-1*reverse(ls)/nLevels; 0; ls/nLevels])
    levels = levels[1:(length(levels) - 1)]
    return levels
end
#
function shubhamBoolean(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int; nLevels::Union{Int, Vector{Int}, String} = 2, 
    vaibhav::Bool = false, turnOffNodes::Array{Int,1} = Int[])
    n_nodes = size(update_matrix,1)
    if (nLevels isa Int)
        nLevels = [nLevels for i in 1:n_nodes]
    elseif (nLevels isa String)
        nLevels = [sum(update_matrix[i,:] .!= 0) + 1 for i in 1:n_nodes]
    end
    sVecList = [getStateVec(nLevels[i]) for i in 1:n_nodes]
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
    updOriginal = copy(update_matrix)  # need it to calculate frustration
    # dividing each column by the corresponding node indegree to ensure that the product is betweenn -1 and 1
    for i in 1:n_nodes
        n = sum(update_matrix[:,i] .!= 0)        
        if n == 0
            n = 1
        end
        update_matrix[:, i] = update_matrix[:, i]/n
    end
    # can't do 2U + I because the previous state has to be retained, which can be any number
    update_matrix2 = sparse(update_matrix')
    stateList = getindex.([rand(sVecList[i], nInit) for i in 1:n_nodes], (1:nInit)')
    for i in 1:nInit
        state = stateList[:,i] #pick random state
        init = stateConvert(state, nLevels)
        flag = 0
        time = 0
        uList = rand(1:n_nodes, nIter)
        states = sVecList
        levels = [getLevels(nLevels[k]) for k in 1:n_nodes]
        j = 1
        while j < nIter
            st = copy(state)
            u = uList[j]
            prod = update_matrix2*st
            s1 = stateChar(prod, state, levels, states, vaibhav)
            if s1[u] != state[u]
                j = j + 1
                time = time + 1
                state[u] = s1[u]
                continue
            end
            while s1[u] == state[u]
                if iszero(j%10) # check after every ten steps,hopefully reduce the time
                    if s1 == state
                        flag = 1
                        break
                    end
                end
                j = j + 1
                time = time + 1
                if j > nIter
                    break
                end
                u = uList[j]
            end
            if flag == 1
                break
            end
            state[u] = s1[u]
        end
        fr = shubhamFrust(state, findnz(sparse(updOriginal)))
        fin = stateConvert(state, nLevels)
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


# function shubhamBoolean(update_matrix::Array{Int,2},
#     nInit::Int, nIter::Int; nLevels::Union{Int, Array{Int, 1}} = 2, 
#     vaibhav::Bool = false, turnOffNodes::Array{Int,1} = Int[])
#     n_nodes = size(update_matrix,1)
#     ls = collect(1:nLevels)
#     sVecList = [getStateVec(nLevels[i]) for i in 1:n_nodes]
#     initVec = []
#     finVec = []
#     flagVec = []
#     frustVec = []
#     timeVec = []
#     if (vaibhav)
#         if (length(turnOffNodes) == 0)
#             vaibhav = [true for i in 1:n_nodes]
#         else
#             vaibhav = [i in turnOffNodes for i in 1:n_nodes]
#         end
#     else
#         vaibhav = [false for i in 1:n_nodes]
#     end
#     # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
#     update_matrix = float(update_matrix)
#     updOriginal = copy(update_matrix)
#     for i in 1:n_nodes
#         n = length([(i,j) for j=1:n_nodes if update_matrix[j,i] == 0])        
#         n = n_nodes - n
#         if n == 0
#             n = 1
#         end
#         update_matrix[:, i] = update_matrix[:, i]/n
#     end
#     update_matrix2 = update_matrix'
#     stateList = getindex.([rand(sVecList[i], nInit) for i in 1:n_nodes], (1:nInit)')
#     @showprogress for i in 1:nInit
#         state = stateList[:,i] #pick random state
#         init = stateConvert(state, nLevels)
#         flag = 0
#         time = 0
#         uList = rand(1:n_nodes, nIter)
#         states = sVecList
#         levels = [getLevels(nLevels[i]) for i in 1:n_nodes]
#         for j in 1:nIter
#             time = time + 1
#             st = copy(state)
#             u = uList[j]
#             st[u] = update_matrix2[u:u,:]*st
#             st[u] = stateChar(st[u], state[u], levels[u], states[u], vaibhav[u])
#             if iszero(j%10 && state[u] == st[u]) # check after every 10 steps,hopefully reduce the time
#                 s2 = stateChar(update_matrix2*st, state, levels, states, vaibhav)
#                 if s2 == state
#                     flag = 1
#                     break
#                 end
#             end
#             state = st
#         end
#         fr = shubhamFrust(state, findnz(sparse(updOriginal)))
#         fin = stateConvert(state, nLevels)
#         push!(frustVec, fr)
#         push!(initVec, init)
#         push!(finVec, fin)
#         push!(flagVec, flag)
#         push!(timeVec, time)     
#         # push!(states_df, (init, fin, flag))
#     end
#     states_df = DataFrame(init=initVec, 
#             fin = finVec, flag = flagVec)
#     frust_df = DataFrame(fin = finVec, 
#         frust = frustVec)
#     frust_df = unique(frust_df, :fin)
#     timeData = DataFrame(fin = finVec, time = timeVec)
#     timeData = groupby(timeData, :fin)
#     timeData = combine(timeData, :time => avg, renamecols = false)
#     frust_df = innerjoin(frust_df, timeData, on = :fin)
#     return states_df, frust_df
# end
