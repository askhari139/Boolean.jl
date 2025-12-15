function checkAttractor()

end

function syncUpdate(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, stateRep::Int, vaibhav::Bool, 
    turnOffNodes::Array{Int,1}, 
    kdNodes::Array{Int, 1}, oeNodes::Array{Int, 1};
    stateList::Vector{Vector{Int}} = Vector{Vector{Int}}()
    )

    n_nodes = size(update_matrix,1)
    density = sum(update_matrix .!= 0) / (n_nodes^2)

    stateVec = stateRep == 0 ? [0, 1] : [-1, 1]

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
    if (length(stateList) == 0)
        stateMat = getindex.([rand(stateVec, nInit) for _ in 1:n_nodes], (1:nInit)')
        if !isempty(kdNodes)
            stateMat[kdNodes, :] .= stateVec[1]
        end
        if !isempty(oeNodes)
            stateMat[oeNodes, :] .= stateVec[2]
        end
        stateList = [stateMat[:, i] for i in 1:nInit]
    end
    nInit = length(stateList)
        # Pre-allocate output vectors

    attractor_list = []
    attractor_size = []
    initVec = Vector{String}(undef, nInit)
    finVec = Vector{String}(undef, nInit)
    flagVec = Vector{Int}(undef, nInit)
    frustVec = Vector{Float64}(undef, nInit)
    timeVec = Vector{Int}(undef, nInit)
    # Precompute nonzero structure for frustration calculation
    nz_structure = findnz(sparse(update_matrix))

    @showprogress for i in 1:nInit
        state = stateList[i]
        init_state = zeroConv(state)
        init = join(init_state, "_")
        flag = 0
        time = 1
        # uList = rand(1:n_nodes, nIter)

        j = 1
        while j <= nIter
            s1 = updFunc(update_matrix2 * state)
            if s1 in attractor_list
                attLength = []
                attractor_list = [attractor_list s1]

            end
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