#=
Author : Kishore Hari
function name : asyncUpdate
Description : simulate a network using ising formalism for multiple random initial conditions
Inputs :
    update_matrix : 2-D array; adjecency matrix of the network. each ij element depicts the edge from node i to node j.
    nInit : integer; number of initial conditions to use for simulation
    nIter : integer; maximum number of time steps for simulation
    
Active outputs :
    state_df : DataFrame; 3 columns - [init, fin, flag]:
        init : string; Initial state
        fin : string; Final state
        flag : integer; 1 - fin is a steady state, 0 - not
    Nodes : string array; List of nodes with the order in which state is listed
Passive outputs :
    The function writes states_df to a csv file if csv = true in input
=#
function asyncUpdate(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int)
    n_nodes = size(update_matrix,1)
    stateVec = Int[-1,1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_matrix2 = sparse(float(update_matrix2'))
    @showprogress for i in 1:nInit
        state = rand(float(stateVec), n_nodes) #pick random state
        init = join(["'", join(Int.(replace(x -> x == -1 ? 0 : x, state))), "'"])
        flag = 0
        time = 0
        for j in 1:nIter
            time = time + 1
            s1 = float(sign.(update_matrix2*state))
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = frustration(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(Int.(replace(x -> x == -1 ? 0 : x, state))), "'"])
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

# function asyncUpdate(update_matrix::Array{Int,2},
#     nInit::Int, nIter::Int)
#     n_nodes = size(update_matrix,1)
#     stateVec = Int[-1,1]
#     initVec = []
#     finVec = []
#     flagVec = []
#     frustVec = []
#     timeVec = []
#     # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
#     update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
#     update_matrix2 = sparse(update_matrix2')
#     @showprogress for i in 1:nInit
#         state = rand(stateVec, n_nodes) #pick random state
#         init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
#         flag = 0
#         time = 0
#         for j in 1:nIter
#             time = time + 1
#             s1 = sign.(update_matrix2*state)
#             u = rand(1:n_nodes, 1)
#             if iszero(j%2) # check after every two steps,hopefully reduce the time
#                 if s1 == state
#                     flag = 1
#                     break
#                 end
#             end
#             state[u] = s1[u]
#         end
#         fr = frustration(state, findnz(sparse(update_matrix)))
#         fin = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
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
#     # print(frust_df)
#     return states_df, frust_df
# end


function asyncUpdate2(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int)
    n_nodes = size(update_matrix,1)
    stateVec = Int[0,1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_matrix2 = sparse(update_matrix2')
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(state), "'"])
        flag = 0
        time = 0
        for j in 1:nIter
            time = time + 1
            s1 = zeroConv(update_matrix2*state)
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = frustration0(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(state), "'"])
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


function async_stg(update_matrix::Array{Int,2})
    nSpecies = size(update_matrix)[1]
    state_list = liststates_df(nSpecies)
    stg = DataFrame(Init=String[], Fin = String[])
    for i in state_list
        i1 = sign.(i'*update_matrix)
        for j in 1:nSpecies
            si = join(["'", join(replace(x -> x == -1 ? 0 : x, i)), "'"])
            sj = copy(i)
            sj[j] = i1[j]
            sj = join(["'", join(replace(x -> x == -1 ? 0 : x, sj)), "'"])
            push!(stg, (si, sj))
        end
    end
    return stg
end

function asyncRandUpdate(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, randVec::Array{Float64,1})
    n_nodes = size(update_matrix,1)
    nzId = enumerate(findall(update_matrix.!=0))
    update_matrix = float(update_matrix)
    for (i,j) in nzId
        update_matrix[j] = update_matrix[j]*randVec[i]
    end
    stateVec = [-1.0,1.0]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    # minVal = minimum([abs(update_matrix[j]) for (i,j) in nzId])
    # update_matrix2 = update_matrix + Matrix(I, n_nodes, n_nodes)*(minVal/2)
    update_matrix2 = sparse(update_matrix')
    for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(Int.(replace(x -> x == -1.0 ? 0 : 1, state))), "'"])
        flag = 0
        time = 0
        for j in 1:nIter
            time = time + 1
            s1 = float(sign.(update_matrix2*state))
            s1 = [s1[i] == 0 ? state[i] : s1[i] for i in 1:n_nodes]
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            if (s1[u] != 0)
                state[u] = s1[u]
            end
        end
        fr = frustration(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(Int.(replace(x -> x == -1.0 ? 0 : 1, state))), "'"])
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


function asyncRandUpdate0(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, randVec::Array{Float64,1})
    n_nodes = size(update_matrix,1)
    nzId = enumerate(findall(update_matrix.!=0))
    update_matrix = [Float64(i) for i in update_matrix]
    for (i,j) in nzId
        update_matrix[j] = update_matrix[j]*randVec[i]
    end
    stateVec = Int[0,1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    minVal = minimum([abs(update_matrix[j]) for (i,j) in nzId])
    update_matrix2 = update_matrix + Matrix(I, n_nodes, n_nodes)*(minVal/2)
    update_matrix2 = sparse(update_matrix2')
    for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(state), "'"])
        flag = 0
        for j in 1:nIter
            time = time + 1
            s1 = zeroConv(update_matrix2*state)
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = frustration(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(state), "'"])
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

function asyncOEDE(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, OEID::Union{Int,Array{Int,1}}=[], 
    DEID::Union{Int, Array{Int,1}}=[])
    n_nodes = size(update_matrix,1)
    if length(OEID) != 0
        update_matrix[:,OEID] = repeat([0], n_nodes, length(OEID))
    end
    if length(DEID) != 0
        update_matrix[:,DEID] = repeat([0], n_nodes, length(OEID))
    end
    stateVec = Int[-1,1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_matrix2 = sparse(update_matrix2')
    for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        if length(OEID) != 0
            state[OEID] = 1
        end
        if length(DEID) != 0
            state[DEID] = -1
        end
        init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        flag = 0
        for j in 1:nIter
            time = time + 1
            s1 = sign.(update_matrix2*state)
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = frustration(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
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

