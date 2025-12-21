# ============================
# Random initial state generator
# ============================
"""
    random_state(N, n; rng=Random.GLOBAL_RNG)

Return a random integer state vector of length `N` with elements in `-n:n` **excluding 0**.
Each element represents a discrete state `k/n`.
"""
function random_state(N::Integer, n::Integer; rng=Random.GLOBAL_RNG, stateRep::Int = -1)
    if stateRep == 0
        return rand(rng, [k for k in -n:n if k != 0], N)
    else
        return rand(rng, [k for k in -n:n if k != 0], N)
    end
end

"""
    encode_state(s, n)

Encode a state vector `s` (integers in -n..-1,1..n, no zeros) into a string.

- Negative states -n..-1 → 0..(n-1)
- Positive states 1..n   → n..(2n-1)
- Result: string of mapped integers joined by "_"
"""
function encode_state(s::Vector{Int}, n::Int)
    mapped = Vector{Int}(undef, length(s))
    @inbounds for i in 1:length(s)
        val = s[i]
        if val < 0
            mapped[i] = val + n   # e.g. -n → 0, -1 → n-1
        else
            mapped[i] = val + n - 1   # 1 → n, n → 2n-1
        end
    end
    return join(mapped, "_")
end


# ============================
# Main asynchronous simulation function
# ============================
"""
    simulate_async!(s, J, n; steps=1000, rng, stateList, record_every,
                    perturb_nodes, perturb_mode, perturb_interval)

Simulate the multilevel Boolean network asynchronously **in-place** on integer-scaled state `s` (`-n..n`, no zeros).

# Arguments
- `s::Vector{Int}` : initial state vector (length = number of nodes), modified in-place
- `J::AbstractMatrix{Int}` : adjacency/sign matrix where `J[i,j]` is the influence of node `i` on node `j`
- `n::Int` : denominator (levels are `k/n`, k=-n..n, no zeros)

# Keyword Arguments
- `steps::Int=1000` : total number of asynchronous update steps
- `rng` : random number generator (default `Random.GLOBAL_RNG`)
- `stateList::Vector{Vector{Int}}=Vector{Vector{Int}}()` : optional history recording
- `record_every::Int=1` : record every `record_every` steps
- `perturb_nodes::Int=0` : number of nodes to perturb at each perturbation event
- `perturb_mode::Symbol=:mirror` : `:mirror` flips `k → -k`, `:full` flips to extreme `-n*sign(k)`
- `perturb_interval::Int=typemax(Int)` : steps between perturbations (`Inf` = no perturbations)

# Returns
`(s, stateList)` : final state vector and recorded history
"""
function simulate_async!(
    s::Vector{Int},
    J::AbstractMatrix{Int},
    n::Int;
    steps::Int = 1000,
    rng = Random.GLOBAL_RNG,
    record_every::Int = 1,
    perturb_nodes::Int = 0,
    perturb_mode::Symbol = :mirror,
    perturb_interval::Int = typemax(Int)
)
    N = length(s)
    @assert size(J,1) == N && size(J,2) == N "J must be square N×N matching s"

    # Precompute denominators: d[i] = sum_j |J[j,i]|
    stateList = Vector{Vector{Int}}()
    push!(stateList, copy(s))
    d = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        ssum = 0
        for j in 1:N
            ssum += abs(J[j,i])
        end
        d[i] = ssum
    end

    # Main asynchronous loop
    @inbounds for t in 1:steps
        # pick one random node
        i = rand(rng, 1:N)
        di = d[i]
        if di > 0
            # compute input sum
            M = 0
            @inbounds for j in 1:N
                M += J[j,i] * s[j]
            end

            if M != 0
                if M < -(n-1)*di
                    newk = -n
                elseif M > (n-1)*di
                    newk = n
                else
                    if M < 0
                        k = (-M + di - 1) ÷ di
                        newk = -k
                    else
                        k = (M + di - 1) ÷ di
                        newk = k
                    end
                end
                # enforce non-zero state
                if newk == 0
                    newk = sign(M)  # ±1
                end
                s[i] = newk
            end
        end

        # Apply perturbation if scheduled
        if (perturb_interval > 0 && (t % perturb_interval == 0) || t == 10) && perturb_nodes > 0
            idxs = rand(rng, 1:N, perturb_nodes)  # sample with replacement
            for i in idxs
                if perturb_mode == :mirror
                    s[i] = -s[i]
                elseif perturb_mode == :full
                    s[i] = -n * sign(s[i])
                elseif perturb_mode == :next
                    s[i] = s[i] + rand([-1,1])
                else
                    error("Unknown perturbation mode: $perturb_mode")
                end
            end
        end

        # Record state
        if (t % record_every == 0)
            push!(stateList, copy(s))
        end
    end
    times = collect(0:(length(stateList)-1)) .* record_every
    return [s ./ n for s in stateList], times
end

# ============================
# Non-mutating wrapper
# ============================
"""
    simulate_async(J, n; s0, steps, rng, stateList, record_every,
                   perturb_nodes, perturb_mode, perturb_interval)

Non-mutating wrapper that copies `s0` or generates a random initial state, then calls `simulate_async!`.
"""
function simulate_async(
    J::AbstractMatrix{Int},
    n::Int;
    s0::Union{Nothing,AbstractVector{Int}} = nothing,
    steps::Int = 1000,
    rng = Random.GLOBAL_RNG,
    record_every::Int = 1,
    perturb_nodes::Int = 0,
    perturb_mode::Symbol = :mirror,
    perturb_interval::Int = typemax(Int)
)
    N = size(J,1)
    @assert size(J,2) == N "J must be square N×N"
    s_init = s0 === nothing ? random_state(N, n; rng=rng) : copy(collect(s0))
    return simulate_async!(
        s_init, J, n;
        steps=steps, rng=rng, record_every=record_every,
        perturb_nodes=perturb_nodes, perturb_mode=perturb_mode, perturb_interval=perturb_interval
    )
end


"""
    simulate_multiple_states_to_df(states, J, names; kwargs...)

Run `simulate_async` on each state in `states` and store all recorded states in a DataFrame.

# Arguments
- `states::Vector{Vector{Int}}` : list of initial states
- `J::AbstractMatrix{Int}` : adjacency/sign matrix
- `names::Vector{String}` : names of nodes, length must match state vector length

# Keyword Arguments
- All named arguments of `simulate_async` can be passed here (steps, perturb_nodes, etc.)

# Returns
- `df::DataFrame` : each row is a recorded state, columns are node names, plus `:run` column indicating which initial state
"""
function simulate_multiple_states_to_df(topoFile::String, n::Int;
                                        states::Vector{Vector{Int}}=Vector{Vector{Int}}(),
                                        nSim::Int=10,
                                        kwargs...)
    J,nms = topo2interaction(topoFile)
    N = length(nms)
    
    all_rows = DataFrame()
    if isempty(states)
        states = [random_state(N, n) for _ in 1:10]
    end
    for (run_idx, s0) in enumerate(states)
        for i in 1:nSim
        # Run simulation
            history, times = simulate_async(J, n; s0=s0, kwargs...)
            # Convert history (Vector of Vectors) to DataFrame
            if !isempty(history)
                df_run = DataFrame(hcat(history...)', Symbol.(nms))
                df_run[!, :Iteration] = fill(i, size(df_run, 1))
                df_run[!,:state_id] = fill(run_idx, size(df_run, 1))  # optional: keep track of which initial state
                df_run[!,:Time] = times
                append!(all_rows, df_run)
            end
        end
    end
    
    return all_rows
end


function simulate_sync_logical(x_init, F, N; record_every::Int = 1,
    perturb_nodes::Int = 0,
    perturb_interval::Int = typemax(Int), max_steps::Int = 1000)
    x = copy(x_init)
    trajectory = Dict{Vector{Int}, Int}()
    trajectory[copy(x)] = 0
    
    found_attractor = false
    cont_pert = perturb_interval > 0 && perturb_interval != Inf && perturb_nodes != 0
    traj = Vector{Vector{Int}}()
    for step in 1:max_steps
        # Synchronous update
        x_new = [F[i](x) for i in 1:N]
        
        # Check if we've closed a cycle in THIS trajectory
        if !cont_pert && haskey(trajectory, x_new)
            cycle_start = trajectory[x_new]
            cycle_states = sort(collect(keys(trajectory)), by=s -> trajectory[s])
            cycle_states = cycle_states[(cycle_start + 1):end]
            if (cycle_states != [x_new])
                push!(cycle_states, x_new)
            end
            # Check if this cycle matches an existing attractor
            found_attractor = true
            break
        end
        
        trajectory[copy(x_new)] = step
        x = x_new
        if cont_pert && (step%perturb_interval == 0 || step == 10)
            # print("enter the dragon")
            if step%record_every == 0
                # print("save the dragon")
                push!(traj, x)
            end
            idxs = rand(1:N, perturb_nodes)  # sample with replacement
            for i in idxs
                if x[i] == 0
                    k = 1
                else
                    k = 0
                end
                x[i] = k
            end
        end
    end
    
    # Max steps reached without finding attractor
    if !found_attractor
        cycle_states = [x]
    end
    if cont_pert
        # println("test")
        return traj
    else
        return cycle_states
    end
end

"""
    simulate_multiple_states_to_df(states, J, names; kwargs...)

Run `simulate_async` on each state in `states` and store all recorded states in a DataFrame.

# Arguments
- `states::Vector{Vector{Int}}` : list of initial states
- `J::AbstractMatrix{Int}` : adjacency/sign matrix
- `names::Vector{String}` : names of nodes, length must match state vector length

# Keyword Arguments
- All named arguments of `simulate_async` can be passed here (steps, perturb_nodes, etc.)

# Returns
- `df::DataFrame` : each row is a recorded state, columns are node names, plus `:run` column indicating which initial state
"""
function simulate_multiple_states_to_df_logical(rules_file::String; 
                                        states::Vector{Vector{Int}}=Vector{Vector{Int}}(),
                                        nSim::Int=10,
                                        kwargs...)
    F, I, N, degrees, variables, constants = getNodeFunctions(rules_file)
    nms = vcat(variables, constants)
    all_rows = DataFrame()
    if isempty(states)
        attrFile = replace(rules_file, ".txt" => "_attractors.csv")
        if isfile(attrFile)
            df = CSV.read(attrFile, DataFrame)
            n_states = min(size(df, 1), 10)
            states = vcat(
                split.(df.states[1:div(n_states, 2)], ";"),
                split.(df.states[div(n_states, 2)+1:end], ";")
            )
            states = [[parse(Int, c) for c in split(state[1], "_")] for state in states]
        else
            states = [random_state(N, 1;stateRep = 0) for _ in 1:10]
        end
    end
    for (run_idx,s0) in enumerate(states)
        for i in 1:nSim
        # Run simulation
            history = simulate_sync_logical(s0, F, N; kwargs...)
            # Convert history (Vector of Vectors) to DataFrame
            if !isempty(history)
                df_run = DataFrame(hcat(history...)', Symbol.(nms))
                df_run[!, :Iteration] = fill(i, size(df_run, 1))
                df_run[!, :state_next] = [format_attractor_states([h]) for h in history]
                df_run[!,:state_init] = fill(format_attractor_states([s0]), size(df_run, 1))  
                df_run[!, :state_id] = fill(run_idx, size(df_run, 1))# optional: keep track of which initial state
                append!(all_rows, df_run)
            end
        end
    end
    
    return all_rows
end