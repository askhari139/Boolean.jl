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

function update_function(F::Vector{Function}, I::Vector{Vector{Int}}, 
    N::Int, mode::Symbol, X::Vector{Int})
    return update_dnf(F, N, X)
end

function update_function(F::Matrix, I::Vector{Vector{Int}}, 
    N::Int, mode::Symbol, X::Vector{Int})
    if mode == :ising
        return update_ising(F,X)
    elseif mode == :nising
        return update_nising(F,X)
    end
end

function update_function(F::Vector{Vector{Int}}, I::Vector{Vector{Int}}, 
    N::Int, mode::Symbol, X::Vector{Int})
    return update_truthtable(F, I, N, X)
end


function get_synchronous_stg(F::Union{Vector{Function}, Vector{Vector{Int}}, Matrix}, 
        I::Vector{Vector{Int}}; 
        exact::Bool = false, mode::Symbol=:logical, nsim::Int = 100000)
    N = size(F,1)
    stg = Dict{Vector{Int64}, Vector{Int64}}()
    if N < 15
        exact = true
    end
    if exact
        initial_conditions = [collect(digits(i, base=2, pad=N)) for i in 0:(2^N - 1)]
    else
        initial_conditions = [rand([0, 1], N) for _ in 1:nsim]
    end
    if mode == :ising
        initial_conditions = [Int.((i.-0.5)./0.5) for i in initial_conditions]
    end
    if !(mode in [:logical, :ising, :nising])
        return stg
    end
    for X in initial_conditions
        y = update_function(F,I, N,mode,X)
        if mode == :ising
            X = Int.((X.+1)./2)
        end
        # print(X)
        # stg[Int.((X.+1)./2)] = Int.((Boolean.update_function(F,I, N, X) .+1)./2)
        stg[X] = y
    end
    return stg
end

function topo_to_stg(topo_file::String; mode = :ising, kwargs...)
    update_matrix, nodes = topo2interaction(topo_file)
    I = Vector{Vector{Int}}()
    N = length(nodes)
    # update_matrix = 2 * update_matrix' + Matrix(I, n_nodes, n_nodes)
    stg = get_synchronous_stg(update_matrix, I; mode = mode, kwargs...)
    return stg
end

function bool_to_stg(bool_file::String; kwargs...)
    F,I,N,_ = getNodeFunctions(bool_file)
    stg = get_synchronous_stg(F, I; mode = :logical, kwargs...)
    return stg
end



function stg_to_digraph(stg::Dict)
    states = collect(keys(stg))
    idx = Dict(s => i for (i, s) in enumerate(states))

    g = DiGraph(length(states))

    for (s, s′) in stg
        add_edge!(g, idx[s], idx[s′])
    end

    return g, states, idx
end

function attractors(g::DiGraph, sccs)
    attractor_sccs = []

    for comp in sccs
        is_terminal = true

        for v in comp
            for w in outneighbors(g, v)
                w ∈ comp || (is_terminal = false; break)
            end
            is_terminal || break
        end

        is_terminal && push!(attractor_sccs, comp)
    end

    attractor_sccs
end



function get_attractors_from_stg(stg)
    g, states, idx = stg_to_digraph(stg)
    sccs = strongly_connected_components(g)
    comp_id = zeros(Int, nv(g))
    for (i, comp) in enumerate(sccs)
        for v in comp
            comp_id[v] = i
        end
    end
    
    # Find terminal SCCs
    is_terminal = trues(length(sccs))
    for (i, comp) in enumerate(sccs)
        for v in comp, w in outneighbors(g, v)
            if comp_id[w] != i
                is_terminal[i] = false
                break
            end
        end
    end
    
    # Compute basins via reverse BFS
    reverse_g = reverse(g)
    results = Dict{Vector{eltype(states)}, Float64}()
    
    for (i, comp) in enumerate(sccs)
        is_terminal[i] || continue
        visited = falses(nv(g))
        queue = copy(comp)
        visited[comp] .= true
        while !isempty(queue)
            v = popfirst!(queue)
            for w in outneighbors(reverse_g, v)
                visited[w] && continue
                visited[w] = true
                push!(queue, w)
            end
        end
        results[states[comp]] = count(visited) / nv(g)
    end
    
    return results
end