##################################################
# Not using this script anymore
##################################################


# function numChar(x, s0, levels, states, vaibhav)
#     if x == 0
#         if vaibhav
#             y = 0
#         else
#             y = s0
#         end
#     else
#         compares = sum(x .>= levels)
#         if compares == 0
#             compares = 1
#         end
#         y = states[compares]
#     end
# end


# function stateChar(state::AbstractArray, s0, levels, states, vaibhav)
#     for i in 1:length(state)
#         x = state[i]
#         y = numChar(x, s0[i], levels[i], states[i], vaibhav)
#         state[i] = y
#     end
#     return state
# end

# function shubhamFrust(state::Array{Float64,1}, 
#     nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}})
#     frustration = 0
#     nEdges = length(nonZeros[1])
#     for (x,y,v) in zip(nonZeros...)
#         s = state[x]*state[y]*v
#         # if s<0
#         #     if state[x] == state[y]
#         #         frustration = frustration + 1
#         #     else
#         #         frustration = frustration + 0.5
#         #     end
#         # end
#         if s < 0
#             frustration = frustration + abs(state[x]*state[y])
#         end
#     end
#     frustration = frustration/nEdges
#     return frustration
# end

# function stateConvert(state, nLevels)
#     state = Int.(state.*nLevels)
#     for i in eachindex(state)
#         if state[i] < 0
#             state[i] = state[i] + nLevels[i]
#         else
#             state[i] = state[i] + nLevels[i] - 1
#         end
#     end
#     state = join(state, "_")
#     # end
#     # state = join(["'", join(replace(x -> x < 0 ? x+nLevels : x+nLevels-1, state)), "'"])
#     return state
# end

# function getStateVec(nLevels::Int = 2)
#     ls = collect(1:nLevels)
#     stateVec = [-1*reverse(ls)/nLevels; ls/nLevels]
#     return stateVec
# end

# function getLevels(nLevels::Int = 2)
#     ls = collect(1:nLevels)
#     levels = [-1*reverse(ls)/nLevels; 0; ls/nLevels]
#     levels = levels[1:(length(levels) - 1)]
#     return levels
# end

# function shubhamBoolean(update_matrix::Array{Int,2},
#     nInit::Int, nIter::Int; nLevels::Union{Int, Array{Int,1}} = 2, 
#     vaibhav::Bool = false)
#     n_nodes = size(update_matrix,1)
#     if (typeof(nLevels) == Int)
#         nLevels = [nLevels for i in 1:n_nodes]
#     end
#     sVecList = [getStateVec(nLevels[i]) for i in 1:n_nodes]
#     initVec = []
#     finVec = []
#     flagVec = []
#     frustVec = []
#     timeVec = []
#     # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
#     update_matrix = float(update_matrix)
#     for i in 1:n_nodes
#         n = length([(i,j) for j=1:n_nodes if update_matrix[j,i] == 0])        
#         n = n_nodes - n
#         if n == 0
#             n = 1
#         end
#         update_matrix[:, i] = update_matrix[:, i]/n
#     end
#     # update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
#     # nzId = enumerate(findall(update_matrix.!=0))
#     # minVal = minimum([abs(update_matrix[j]) for (i,j) in nzId])
#     # update_matrix2 = update_matrix + Matrix(I, n_nodes, n_nodes)*(minVal/2)
#     update_matrix2 = sparse(update_matrix')
#     stateList = getindex.([rand(sVecList[i], nInit) for i in 1:n_nodes], (1:nInit)')
#     @showprogress for i in 1:nInit
#         state = stateList[:,i] #pick random state
#         init = stateConvert(state, nLevels)
#         flag = 0
#         time = 0
#         levels = [getLevels(nLevels[i]) for i in 1:n_nodes]
#         for j in 1:nIter
#             time = time + 1
#             st = copy(state)
#             s1 = stateChar(update_matrix2*st, state, levels, sVecList, vaibhav)
#             if iszero(j%2) # check after every two steps,hopefully reduce the time
#                 if s1 == state
#                     flag = 1
#                     break
#                 end
#             end
#             #update a random node
#             u = rand(1:n_nodes, 1)
#             state[u] = s1[u]
#         end
#         fr = shubhamFrust(state, findnz(sparse(update_matrix)))
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
