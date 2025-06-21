"""
    asyncUpdate(update_matrix::Array{Int,2},
                nInit::Int, nIter::Int, stateRep::Int, vaibhav::Bool, 
                turnOffNodes::Array{Int,1}, 
                kdNodes::Array{Int,1}, oeNodes::Array{Int,1})

Perform asynchronous updates on a network based on the provided update matrix and settings.

# Arguments
- `update_matrix::Array{Int,2}`: The adjacency matrix of the network (values: -1, 0, 1).
- `nInit::Int`: Number of initial random states to simulate.
- `nIter::Int`: Maximum number of asynchronous updates per simulation.
- `stateRep::Int`: Type of state representation (`0` for {0,1} states, else {-1,1} states).
- `vaibhav::Bool`: Whether to "turn off" nodes dynamically (special rule handling).
- `turnOffNodes::Array{Int,1}`: List of nodes to apply the Vaibhav turn-off rule.
- `kdNodes::Array{Int,1}`: List of nodes forced into constant knockdown (-1 or 0 depending on `stateRep`).
- `oeNodes::Array{Int,1}`: List of nodes forced into constant overexpression (+1).

# Returns
- `states_df::DataFrame`: A table containing the initial and final network states (`init`, `fin`) and a flag (`flag`) indicating convergence.
- `frust_df::DataFrame`: A table containing each unique final state (`fin`), its network frustration score (`frust`), and the average time steps (`time`) taken to reach that state.

# Notes
- Asynchronous update: At each step, a random node is selected and updated.
- Convergence is checked every 10 steps to reduce computational load.
- Forced nodes (`kdNodes`, `oeNodes`) are never updated.
- Frustration quantifies how "unsatisfied" the network is (lower is better).
- Uses different conventions depending on `stateRep`:
    - If `stateRep == 0`, a (0,1) binary system is used.
    - Otherwise, a (-1,1) spin-like system is used.
- Supports sparse optimization if the network is large and sparse.

# Example
```julia
states_df, frust_df = asyncUpdate(update_mat, 100, 1000, 1, false, Int[], Int[], Int[])
"""
function asyncUpdate(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, stateRep::Int, vaibhav::Bool, 
    turnOffNodes::Array{Int,1}, 
    kdNodes::Array{Int, 1}, oeNodes::Array{Int, 1})

    n_nodes = size(update_matrix,1)
    density = sum(update_matrix .!= 0) / (n_nodes^2)

    stateVec = stateRep == 0 ? [0, 1] : [-1, 1]

    # Pre-allocate output vectors
    initVec = Vector{String}(undef, nInit)
    finVec = Vector{String}(undef, nInit)
    flagVec = Vector{Int}(undef, nInit)
    frustVec = Vector{Float64}(undef, nInit)
    timeVec = Vector{Int}(undef, nInit)

    # Handle "vaibhav" adjustment
    idMat = Matrix(I, n_nodes, n_nodes)
    if vaibhav
        if isempty(turnOffNodes)
            turnOffNodes = collect(1:n_nodes)
        end
        idMat[turnOffNodes, :] .= 0
    end

    # Preprocess update matrix
    update_matrix2 = 2 * update_matrix + idMat
    if n_nodes > 500 && density < 0.1
        update_matrix2 = sparse(update_matrix2')
    else
        update_matrix2 = update_matrix2'
    end

    updFunc = stateRep == 0 ? zeroConv : signVec

    # Initialize random states
    stateMat = getindex.([rand(stateVec, nInit) for _ in 1:n_nodes], (1:nInit)')
    if !isempty(kdNodes)
        stateMat[kdNodes, :] .= stateVec[1]
    end
    if !isempty(oeNodes)
        stateMat[oeNodes, :] .= stateVec[2]
    end
    stateList = [stateMat[:, i] for i in 1:nInit]

    # Precompute nonzero structure for frustration calculation
    nz_structure = findnz(sparse(update_matrix))

    @showprogress for i in 1:nInit
        state = stateList[i]
        init_state = zeroConv(state)
        init = join(init_state, "_")

        flag = 0
        time = 1
        uList = rand(1:n_nodes, nIter)

        j = 1
        while j <= nIter
            s1 = updFunc(update_matrix2 * state)
            u = uList[j]

            # Skip update if node is knocked down or overexpressed
            if u in kdNodes || u in oeNodes
                j += 1
                time += 1
                continue
            end

            # Perform update if necessary
            if s1[u] != state[u]
                state[u] = s1[u]
                j += 1
                time += 1
                continue
            end

            # Otherwise, wait for change or convergence
            while s1[u] == state[u]
                if iszero(j % 10) && s1 == state
                    flag = 1
                    break
                end
                j += 1
                time += 1
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

        # Calculate frustration
        fr = stateRep == 0 ? frustration(state, nz_structure; negConv = true) : frustration(state, nz_structure)

        fin_state = zeroConv(state)
        fin = join(fin_state, "_")

        initVec[i] = init
        finVec[i] = fin
        flagVec[i] = flag
        frustVec[i] = fr
        timeVec[i] = time
    end

    # Prepare output dataframes
    states_df = DataFrame(init=initVec, fin=finVec, flag=flagVec)
    frust_df = DataFrame(fin=finVec, frust=frustVec)
    frust_df = unique(frust_df, :fin)

    timeData = DataFrame(fin=finVec, time=timeVec)
    timeData = DataFrames.groupby(timeData, :fin)
    timeData = DataFrames.combine(timeData, :time => avg, renamecols=false)

    frust_df = DataFrames.innerjoin(frust_df, timeData, on=:fin)

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
            state_df, frust_df = asyncUpdate(update_matrix, nInit, nIter, stateRep,false,[],[],[])
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
    timeData = DataFrames.groupby(timeData, :fin)
    timeData = DataFrames.combine(timeData, :time => avg, renamecols = false)
    frust_df = DataFrames.innerjoin(frust_df, timeData, on = :fin)
    return states_df, frust_df
end


function asyncRandCont(update_matrix::Union{Array{Int,2}, Array{Float64,2}},
    nInit::Int, nIter::Int, stateRep::Int; randVec::Array{Float64,1} = [0.0], 
    weightFunc::Function = defaultWeightsFunction(0.01), 
    frequency::Int = 1, steadyStates::Bool = true)
    n_nodes = size(update_matrix,1)
    if steadyStates
        updm = sign.(update_matrix)
        states_df, frust_df = asyncUpdate(updm, 10000, 1000, stateRep,false,[],[],[])
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


