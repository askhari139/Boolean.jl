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
    nInit::Int, nIter::Int, stateRep::Int, vaibhav::Bool, 
    turnOffNodes::Array{Int,1})
    n_nodes = size(update_matrix,1)
    stateVec = ifelse(stateRep == 0, [0,1], [-1,1])
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    idMat = Matrix(I, n_nodes, n_nodes)
    if vaibhav
        if length(turnOffNodes) == 0
            turnOffNodes = 1:n_nodes
        end
        idMat[turnOffNodes, :] .= 0
    end

    update_matrix2 = 2*update_matrix + idMat
    update_matrix2 = sparse(update_matrix2')
    updFunc = ifelse(stateRep == 0, zeroConv, signVec)
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(zeroConv(state), "_")
        flag = 0
        time = 1
        uList = rand(1:n_nodes, nIter)
        j = 1
        while j <= nIter
            s1 = updFunc(update_matrix2*state)
            u = uList[j]
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
        if stateRep == 0
            fr = frustration(state, findnz(sparse(update_matrix)); negConv = true)
        else
            fr = frustration(state, findnz(sparse(update_matrix)))
        end
        fin = join(zeroConv(state), "_")
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

function asyncRandUpdate(update_matrix::Union{Array{Int,2}, Array{Float64,2}},
    nInit::Int, nIter::Int, randVec::Array{Float64,1}, stateRep::Int64)
    n_nodes = size(update_matrix,1)
    if typeof(update_matrix) == Array{Int, 2}
        nzId = enumerate(findall(update_matrix.!=0))
        if randVec == [0.0]
            state_df, frust_df = asyncUpdate(update_matrix, nInit, nIter, stateRep)
            return state_df, frust_df
        end
        update_matrix = float(update_matrix)
        for (i,j) in nzId
            update_matrix[j] = update_matrix[j]*randVec[i]
        end
    end
    stateVec = ifelse(stateRep == 0, [0.0,1.0], [-1.0,1.0])
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
        init = join(Int.(zeroConv(state)), "_")
        flag = 0
        time = 1
        uList = rand(1:n_nodes, nIter)
        updFunc = ifelse(stateRep == 0, zeroConv, signVec)
        for j in 1:nIter
            s1 = float(updFunc(update_matrix2*state))
            s1 = [s1[i] == 0 ? state[i] : s1[i] for i in 1:n_nodes]
            u = uList[j]
            while s1[u] == state[u]
                if iszero(j%10) # check after every two steps,hopefully reduce the time
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
        fr = ifelse(stateRep == 0, frustration(state, findnz(sparse(update_matrix)); negConv = true), frustration(state, findnz(sparse(update_matrix)); negConv = false))
        fin = join(Int.(zeroConv(state)), "_")
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


function asyncRandCont(update_matrix::Union{Array{Int,2}, Array{Float64,2}},
    nInit::Int, nIter::Int, stateRep::Int; randVec::Array{Float64,1} = [0.0], 
    weightFunc::Function = defaultWeightsFunction(0.01), 
    frequency::Int = 1, steadyStates::Bool = true)
    n_nodes = size(update_matrix,1)
    if steadyStates
        updm = sign.(update_matrix)
        states_df, frust_df = asyncUpdate(updm, 10000, 1000, stateRep)
        states = states_df[:, :fin]
        states = rand(states, nInit)
    end
    update_matrix = update_matrix'
    nzId = enumerate(findall(update_matrix.!=0))
    if typeof(update_matrix) == Adjoint{Int64, Matrix{Int64}}
        if randVec == [0.0]
            randVec = rand(length(nzId))
        end
        update_matrix = float(update_matrix)
        for (i,j) in nzId
            update_matrix[j] = update_matrix[j]*randVec[i]
        end
    else
        randVec = [update_matrix[j] for (i,j) in nzId]
    end
    stateVec = ifelse(stateRep == 0, [0.0,1.0], [-1.0,1.0])
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    sListUnique = []
    sListKeys = []
    # create a state matrix of size nInit x nIter
    stateMatrix = zeros(Int, nInit, nIter)
    @showprogress for i in 1:nInit
        # print(i)
        update_matrix2 = update_matrix
        if steadyStates
            state = parse.(Float64, split(states[i], "_"))
        else
            state = rand(stateVec, n_nodes) #pick random state
        end
        init = join(Int.(zeroConv(state)), "_")
        flag = 0
        time = 1
        uList = rand(1:n_nodes, nIter)
        updFunc = ifelse(stateRep == 0, zeroConv, signVec)
        sList = [init]
        for j in 2:nIter
            s1 = float(updFunc(update_matrix2*state))
            if iszero(j%frequency)
                randVec = weightFunc(randVec)
                for (k,l) in nzId
                    update_matrix2[l] = update_matrix[l]*randVec[k]
                end
            end
            s1 = [s1[i] == 0 ? state[i] : s1[i] for i in 1:n_nodes]
            u = uList[j]
            state[u] = s1[u]
            st = join(Int.(zeroConv(state)), "_")
            push!(sList, st)
        end
        sUnique = unique(sList)
        for st in sUnique
            if st in sListKeys
                continue
            end
            # push!(sListKeys, st)
            push!(sListUnique, st)
        end
        # sList is the ID of each state in sListUnique
        sInt = [findfirst(x -> x == st, sListUnique) for st in sList]
        stateMatrix[i,:] = sInt
    end
    return stateMatrix, sListUnique
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

