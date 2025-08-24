"""
    getStateVec(nLevels::Int = 2, kd::Bool = false, oe::Bool = false)

Returns a state vector depending on whether the node is knocked down (kd), overexpressed (oe), or neither.  
- If `kd` is true, returns `[-1.0]`.
- If `oe` is true, returns `[1.0]`.
- Otherwise, returns a symmetric vector of values from `-1` to `1` scaled by `nLevels`.
"""
function getStateVec(nLevels::Int = 2, kd::Bool = false, oe::Bool = false)
    if kd
        return [-1.0]
    elseif oe
        return [1.0]
    else
        ls = collect(1:nLevels)
        stateVec = float([-1*reverse(ls)/nLevels; ls/nLevels])
        return stateVec
    end
end

"""
    getLevels(nLevels::Int = 2)

Generates threshold levels for multilevel state discretization.  
Returns an array of levels from `-1` to `1`, excluding the maximum positive value.
"""
function getLevels(nLevels::Int = 2)
    ls = collect(1:nLevels)
    levels = float([-1*reverse(ls)/nLevels; 0; ls/nLevels])
    levels = levels[1:(length(levels) - 1)]
    return levels
end

"""
    numChar(xnC, s0, levels, states, vaibhav)

Given a continuous node value `xnC`, maps it to a discrete value based on `levels` and `states`.
- If `states` has only one element, returns `s0`.
- If `xnC == 0`, returns `0.0` if `vaibhav` is true, else returns `s0`.
- Otherwise, computes the corresponding discrete value according to how `xnC` compares to `levels`.
"""
function numChar(xnC, s0, levels, states, vaibhav)
    if length(states) == 1
        return s0
    end
    if xnC == 0
        return vaibhav ? 0.0 : s0
    end
    sg = sign(xnC)
    xnC = abs(xnC)
    how_big_is_x = count(l -> xnC > l, levels)
    return sg * states[how_big_is_x]
end

"""
    stateChar!(state::Vector{Float64}, s0, levels, states, vaibhav)

In-place transformation of a state vector using `numChar`.  
Updates each element of `state` based on the corresponding `s0`, `levels`, `states`, and `vaibhav`.
"""
function stateChar!(state::Vector{Float64}, s0, levels, states, vaibhav)
    @inbounds for i in eachindex(state)
        state[i] = numChar(state[i], s0[i], levels[i], states[i], vaibhav[i])
    end
    return state
end

"""
    stateChar!(state::Float64, s0, levels, states, vaibhav)

Applies `numChar` to a single Float64 value and returns the transformed result.
"""
function stateChar!(state::Float64, s0, levels, states, vaibhav)
    return numChar(state, s0, levels, states, vaibhav)
end

"""
    shubhamFrust(state::Vector{Float64}, nonZeros::Tuple{Vector{Int64}, Vector{Int64}, Vector{Float64}})

Computes the frustration of a state based on nonzero edges.  
Frustration is defined as the normalized sum of magnitudes where the product of connected node states and edge weight is negative.
"""
function shubhamFrust(state::Array{Float64,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}})
    frustration = 0.0
    nEdges = length(nonZeros[1])
    @inbounds for (x, y, v) in zip(nonZeros...)
        s = state[x] * state[y] * v
        if s < 0
            frustration += abs(s)
        end
    end
    frustration /= nEdges
    return frustration
end

"""
    stateConvert(state::Vector{Float64}, nLevels::Int)

Converts a continuous state vector to a discrete string representation using a single number of levels.  
- Negative values are shifted up by `nLevels`.
- Positive values are shifted by `nLevels - 1`.
- Zero values are mapped to `2 * nLevels`.
Returns a string like `'2_3_1'`.
"""
function stateConvert(state::AbstractVector, nLevels::Int)
    if (nLevels == 0)
        nLevels = 1
    end
    scaled = Int.(nLevels .* state)
    converted = replace(x -> x < 0 ? x + nLevels : x + nLevels - 1, scaled)
    converted = [scaled[i] == 0 ? 2 * nLevels : converted[i] for i in eachindex(scaled)]
    return join(converted, "_")
end


"""
    stateConvert(state::Vector{Float64}, nLevels::Vector{Int})

Converts a continuous state vector to a discrete string representation, where each node can have a different number of levels.  
Handles negative, positive, and zero values independently per node.
Returns a string like `2_3_1`.
"""
function stateConvert(state::AbstractVector, nLevels::AbstractVector{Int})
    scaled = Int.(state .* nLevels)
    @inbounds for i in eachindex(scaled)
        if scaled[i] < 0
            scaled[i] += nLevels[i]
        elseif scaled[i] > 0
            scaled[i] += nLevels[i] - 1
        else
            scaled[i] = 2 * nLevels[i]
        end
    end
    return join(scaled, "_")
end


"""
    shubhamBoolean(update_matrix::Array{Int,2}, nInit::Int, nIter::Int, nLevels::Union{Int, Vector{Int}, String},
                   vaibhav::Bool, turnOffNodes::Array{Int,1}, kdNodes::Array{Int,1}, oeNodes::Array{Int,1})

Simulates asynchronous Boolean (or multi-level) network updates over multiple initializations.

# Arguments
- `update_matrix::Array{Int,2}`: The adjacency matrix representing the regulatory network. Nonzero entries indicate interactions.
- `nInit::Int`: Number of initial conditions (random initializations) to simulate.
- `nIter::Int`: Maximum number of update steps to run for each initial state.
- `nLevels::Union{Int, Vector{Int}, String}`: Number of discrete levels each node can take. 
  - If `Int`, all nodes have the same number of levels.
  - If `Vector`, node-specific levels.
  - If `"inDegree"`, levels are set based on indegree (number of incoming edges + 1).
- `vaibhav::Bool`: If `true`, selectively "turns off" updates at certain nodes listed in `turnOffNodes`.
- `turnOffNodes::Array{Int,1}`: Indices of nodes to be "turned off" during updates if `vaibhav` is true.
- `kdNodes::Array{Int,1}`: Indices of knockdown (KD) nodes, which are forced to minimum state (-1.0).
- `oeNodes::Array{Int,1}`: Indices of overexpression (OE) nodes, which are forced to maximum state (1.0).

# Returns
- `states_df::DataFrame`: DataFrame containing initial (`init`), final (`fin`) states, and a convergence `flag` (1 if converged early, 0 if not).
- `frust_df::DataFrame`: DataFrame containing final states (`fin`), their frustration score (`frust`), and average convergence time (`time`).

# Notes
- If the network is large and sparse, the function automatically uses a sparse matrix format for updates.
- Frustration quantifies how many edges are "unsatisfied" (i.e., product of connected nodes and edge sign is negative).
- KD and OE nodes retain their fixed states throughout the simulation.
- Random asynchronous updates are used: at each step, a random node is selected and updated based on network inputs.
"""
function shubhamBoolean(
    update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, 
    nLevels::Union{Int, Vector{Int}, String}, 
    vaibhav::Bool, 
    turnOffNodes::Array{Int,1},
    kdNodes::Array{Int,1}, oeNodes::Array{Int,1};
    stateList::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
)
    ### --- Preprocessing input arguments ---
    n_nodes = size(update_matrix, 1)
    density = sum(update_matrix .!= 0) / n_nodes^2

    # Process nLevels
    if nLevels isa Int
        nLevels = fill(nLevels, n_nodes)
    elseif nLevels isa String
        nLevels = [sum(update_matrix[i,:] .!= 0) + 1 for i in 1:n_nodes]
    end

    # Process kdNodes, oeNodes, vaibhavNodes
    kdNodesList = [i in kdNodes for i in 1:n_nodes]
    oeNodesList = [i in oeNodes for i in 1:n_nodes]
    vaibhavNodes = vaibhav ? [i in turnOffNodes for i in 1:n_nodes] : falses(n_nodes)

    # Precompute state vectors and levels
    sVecList = [getStateVec(nLevels[i], kdNodesList[i], oeNodesList[i]) for i in 1:n_nodes]
    levels = [getLevels(nLevels[i]) for i in 1:n_nodes]

    # Prepare update matrix
    updOriginal = copy(update_matrix)
    update_matrix = float.(update_matrix)
    nzUpd = findnz(sparse(update_matrix))
    @inbounds for i in 1:n_nodes
        degree = max(sum(update_matrix[:,i] .!= 0), 1)
        update_matrix[:,i] ./= degree
    end
    update_matrix2 = (n_nodes > 500 && density < 0.1) ? sparse(update_matrix') : update_matrix'

    ### --- Simulation Loop ---
    # Initialize random initial states
    if (length(stateList) == 0)
        stateList = getindex.([rand(sVecList[i], nInit) for i in 1:n_nodes], (1:nInit)')
        stateList = [stateList[:, i] for i in 1:nInit]
    end
    nInit = length(stateList)

        # Prepare outputs
    initVec = Vector{String}(undef, nInit)
    finVec = Vector{String}(undef, nInit)
    flagVec = Vector{Int}(undef, nInit)
    frustVec = Vector{Float64}(undef, nInit)
    timeVec = Vector{Int}(undef, nInit)

    @inbounds for i in 1:nInit
        state = stateList[i]
        init = stateConvert(state, nLevels)
        flag = 0
        time = 0
        uList = rand(1:n_nodes, nIter)

        j = 1
        while j < nIter
            u = uList[j]

            # Skip KD/OE nodes
            if kdNodesList[u] || oeNodesList[u]
                j += 1
                time += 1
                continue
            end

            s1 = update_matrix2 * state
            s1 = stateChar!(s1, state, levels, sVecList, vaibhavNodes)

            # State change
            if s1[u] != state[u]
                state[u] = s1[u]
                j += 1
                time += 1
                continue
            end

            # Check convergence
            while s1[u] == state[u]
                if iszero(j % 10)
                    if s1 == state
                        flag = 1
                        break
                    end
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

        # Save output
        fr = shubhamFrust(state, nzUpd)
        fin = stateConvert(state, nLevels)

        initVec[i] = init
        finVec[i] = fin
        flagVec[i] = flag
        frustVec[i] = fr
        timeVec[i] = time
    end

    ### --- Postprocessing ---
    states_df = DataFrame(init = initVec, fin = finVec, flag = flagVec)
    frust_df = DataFrame(fin = finVec, frust = frustVec)
    frust_df = unique(frust_df, :fin)

    timeData = DataFrame(fin = finVec, time = timeVec)
    timeData = groupby(timeData, :fin)
    timeData = combine(timeData, :time => avg, renamecols = false)

    frust_df = innerjoin(frust_df, timeData, on = :fin)

    return states_df, frust_df
end

