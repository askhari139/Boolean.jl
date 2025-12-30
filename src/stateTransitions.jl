
"""
Calculate which nodes don't satisfy target condition based on node-level constraints.

Arguments:
- s_norm: Normalized state vector (values in [-1,1])
- node_conditions: Vector of functions, one per node, each taking a single value and returning Bool
                   Use `nothing` for nodes without constraints

Returns:
- Vector of node indices that don't satisfy their conditions
"""
function find_violating_nodes(s_norm::Vector{Float64}, node_conditions::Vector{Union{Function, Nothing}})
    violating = Int[]
    for i in 1:length(s_norm)
        if node_conditions[i] !== nothing && !node_conditions[i](s_norm[i])
            push!(violating, i)
        end
    end
    return violating
end


"""
Calculate transition probability from starting state to target condition.

Arguments:
- s_start: Starting state vector
- target_condition: Function that takes normalized state and returns Bool (for final check)
- node_conditions: Vector of per-node condition functions (to identify nodes to perturb)
                   Use `nothing` for unconstrained nodes
- J: Interaction matrix
- n: Number of levels
- perturb_fraction: Fraction of violating nodes to perturb (0 to 1)
- perturb_nodes: Optional vector of specific node indices (overrides automatic detection)
- perturb_mode: :mirror, :full, or :next
- perturb_time: Time step when perturbation occurs
- steps: Maximum simulation steps
- steady_state_window: Steps to check for convergence
- rng: Random number generator

Returns:
- 1.0 if final state satisfies target_condition, 0.0 otherwise
"""
function simulate_transition_to_target(
    s_start::Vector{Int},
    target_condition::Function,
    node_conditions::Vector{Union{Function, Nothing}},
    J::AbstractMatrix{Int},
    n::Int;
    perturb_fraction::Float64 = 1.0,
    perturb_nodes::Union{Nothing, Vector{Int}} = nothing,
    perturb_mode::Symbol = :mirror,
    perturb_time::Int = 10,
    steps::Int = 1000,
    steady_state_window::Int = 50,
    rng = Random.GLOBAL_RNG
)
    N = length(s_start)
    @assert size(J, 1) == N && size(J, 2) == N "J must be square N×N"
    @assert length(node_conditions) == N "node_conditions must have length N"
    
    # Copy starting state
    s = copy(s_start)
    
    # Determine which nodes to perturb
    if perturb_nodes === nothing
        # Find nodes that don't satisfy target condition
        s_normalized = s ./ n
        violating_nodes = find_violating_nodes(s_normalized, node_conditions)
        
        if isempty(violating_nodes)
            # Already satisfies condition, no perturbation needed
            nodes_to_perturb = Int[]
        else
            # Perturb a fraction of the violating nodes
            n_perturb = max(1, round(Int, length(violating_nodes) * perturb_fraction))
            nodes_to_perturb = rand(rng, violating_nodes, n_perturb)
        end
    else
        nodes_to_perturb = perturb_nodes
    end
    
    # Precompute denominators: d[i] = sum_j |J[j,i]|
    d = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        ssum = 0
        for j in 1:N
            ssum += abs(J[j, i])
        end
        d[i] = ssum
    end
    
    # Track recent states for steady state detection
    recent_states = Vector{Vector{Int}}()
    
    # Main simulation loop
    @inbounds for t in 1:steps
        # Apply perturbation at specified time
        if t == perturb_time && !isempty(nodes_to_perturb)
            for i in nodes_to_perturb
                if perturb_mode == :mirror
                    s[i] = -s[i]
                elseif perturb_mode == :full
                    s[i] = -n * sign(s[i])
                elseif perturb_mode == :next
                    curr_val = s[i]
                    step = rand(rng, [-1, 1])
                    new_val = curr_val + step
                    # Keep within bounds [-n, n] and avoid 0
                    new_val = clamp(new_val, -n, n)
                    if new_val == 0
                        new_val = step
                    end
                    s[i] = new_val
                end
            end
        end
        
        # Asynchronous update: pick one random node
        i = rand(rng, 1:N)
        di = d[i]
        
        if di > 0
            # Compute input sum
            M = 0
            for j in 1:N
                M += J[j, i] * s[j]
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
                
                # Enforce non-zero state
                if newk == 0
                    newk = sign(M)
                end
                
                s[i] = newk
            end
        end
        
        # Check for steady state every 10 steps
        if t % 10 == 0
            push!(recent_states, copy(s))
            if length(recent_states) > steady_state_window ÷ 10
                popfirst!(recent_states)
            end
            
            # If all recent states are identical, we've reached steady state
            if length(recent_states) >= steady_state_window ÷ 10
                all_same = true
                first_state = recent_states[1]
                for i in 2:length(recent_states)
                    if recent_states[i] != first_state
                        all_same = false
                        break
                    end
                end
                if all_same
                    break  # Converged to steady state
                end
            end
        end
    end
    
    # Normalize to [-1, 1] and check target condition
    s_normalized = s ./ n
    return target_condition(s_normalized) ? 1.0 : 0.0
end


"""
Calculate empirical transition probability by running multiple simulations.
"""
function calculate_transition_probability(
    s_start::Vector{Int},
    target_condition::Function,
    node_conditions::Vector{Union{Function, Nothing}},
    J::AbstractMatrix{Int},
    n::Int;
    n_trials::Int = 100,
    kwargs...
)
    successes = 0
    for trial in 1:n_trials
        result = simulate_transition_to_target(
            s_start, target_condition, node_conditions, J, n; kwargs...
        )
        successes += Int(result)
    end
    return successes / n_trials
end