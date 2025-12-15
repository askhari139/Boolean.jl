using Random
using DataFrames

"""
    find_attractors_synchronous(
        F::Vector{Function},
        N::Int;
        nsim::Int = 500,
        exact::Bool = false,
        initial_conditions::Vector{Vector{Int}} = Vector{Vector{Int}}(),
        max_steps::Int = 1000,
        seed::Int = -1,
        debug::Bool = false
    ) -> Tuple{DataFrame, Dict{Vector{Int}, Int}}

Find attractors using synchronous updates (deterministic).
Uses state memoization since each state always leads to the same attractor.
"""
function find_attractors_synchronous(
    F::Vector{Function},
    N::Int;
    nsim::Int = 100000,
    exact::Bool = false,
    initial_conditions::Vector{Vector{Int}} = Vector{Vector{Int}}(),
    max_steps::Int = 1000,
    seed::Int = -1,
    debug::Bool = false
)
    # Set random seed
    if seed == -1
        seed = rand(UInt32)
    end
    Random.seed!(seed)
    
    # Generate initial conditions
    if exact
        initial_conditions = [collect(digits(i, base=2, pad=N)) for i in 0:(2^N - 1)]
    elseif isempty(initial_conditions)
        initial_conditions = [rand([0, 1], N) for _ in 1:nsim]
    end
    
    # Data structures
    attractors = Vector{Vector{Vector{Int}}}()
    attractor_times = Vector{Vector{Int}}()
    attractor_dict = Dict{Vector{Int}, Int}()  # State memoization for efficiency
    converged_flags = Int[]
    
    # Process each initial condition
    for (iter, x_init) in enumerate(initial_conditions)
        x = copy(x_init)
        trajectory = Dict{Vector{Int}, Int}()
        trajectory[copy(x)] = 0
        
        found_attractor = false
        
        for step in 1:max_steps
            # Check if we've already seen this state in a previous trajectory
            if haskey(attractor_dict, x)
                attractor_id = attractor_dict[x]
                push!(attractor_times[attractor_id], step)
                
                # Memoize all states in current trajectory
                for (state, _) in trajectory
                    if !haskey(attractor_dict, state)
                        attractor_dict[state] = attractor_id
                    end
                end
                
                found_attractor = true
                break
            end
            
            # Synchronous update
            x_new = [F[i](x) for i in 1:N]
            
            # Check if we've closed a cycle in THIS trajectory
            if haskey(trajectory, x_new)
                cycle_start = trajectory[x_new]
                cycle_states = sort(collect(keys(trajectory)), by=s -> trajectory[s])
                cycle_states = cycle_states[(cycle_start + 1):end]
                if (cycle_states != [x_new])
                    push!(cycle_states, x_new)
                end
                # Check if this cycle matches an existing attractor
                attractor_id = find_equivalent_attractor(cycle_states, attractors)
                
                if isnothing(attractor_id)
                    # New attractor
                    attractor_id = length(attractors) + 1
                    push!(attractors, cycle_states)
                    push!(attractor_times, [step])
                    push!(converged_flags, 1)
                    
                    # Memoize all cycle states
                    for state in cycle_states
                        attractor_dict[state] = attractor_id
                    end
                else
                    push!(attractor_times[attractor_id], step)
                end
                
                # Memoize trajectory states
                for (state, _) in trajectory
                    if !haskey(attractor_dict, state)
                        attractor_dict[state] = attractor_id
                    end
                end
                
                found_attractor = true
                break
            end
            
            trajectory[copy(x_new)] = step
            x = x_new
            
            if debug && step % 100 == 0
                println("Iteration $iter, Step $step")
            end
        end
        
        # Max steps reached without finding attractor
        if !found_attractor
            attractor_id = length(attractors) + 1
            push!(attractors, [x])
            push!(attractor_times, [max_steps])
            push!(converged_flags, 0)
            
            for (state, _) in trajectory
                attractor_dict[state] = attractor_id
            end
        end
    end
    
    # Calculate statistics
    basin_sizes = [count(==(id), values(attractor_dict)) for id in 1:length(attractors)]
    basin_proportions = basin_sizes ./ length(initial_conditions)
    avg_times = [sum(times) / length(times) for times in attractor_times]
    
    
    df = DataFrame(
        states = attractors,
        basin_size = basin_proportions,
        time = avg_times,
        flag = converged_flags
    )
    
    return df, attractor_dict
end

"""
    find_attractors_asynchronous(
        F::Vector{Function},
        N::Int;
        nsim::Int = 500,
        initial_conditions::Vector{Vector{Int}} = Vector{Vector{Int}}(),
        max_steps::Int = 1000,
        seed::Int = -1,
        debug::Bool = false
    ) -> DataFrame

Find attractors using asynchronous updates (stochastic).
Each trajectory is independent since dynamics are non-deterministic.
Basin sizes represent the proportion of trajectories reaching each attractor.
"""
function find_attractors_asynchronous(
    F::Vector{Function},
    N::Int;
    nsim::Int = 100000,
    initial_conditions::Vector{Vector{Int}} = Vector{Vector{Int}}(),
    max_steps::Int = 1000,
    seed::Int = -1,
    debug::Bool = false
)
    # Set random seed
    if seed == -1
        seed = rand(UInt32)
    end
    Random.seed!(seed)
    
    # Generate initial conditions
    if isempty(initial_conditions)
        initial_conditions = [rand([0, 1], N) for _ in 1:nsim]
    end
    
    # Data structures - NO state dictionary for async
    attractors = Vector{Vector{Vector{Int}}}()
    attractor_counts = Int[]  # Count trajectories reaching each attractor
    attractor_times = Vector{Vector{Int}}()
    converged_flags = Int[]
    
    # Process each trajectory independently
    for (iter, x_init) in enumerate(initial_conditions)
        x = copy(x_init)
        trajectory = Dict{Vector{Int}, Int}()
        trajectory[copy(x)] = 0
        
        found_attractor = false
        
        for step in 1:max_steps
            # Asynchronous update: pick random node that wants to change
            x_new = copy(x)
            nodes_to_update = randperm(N)
            for i in nodes_to_update
                new_val = F[i](x_new)
                if x_new[i] != new_val
                    x_new[i] = new_val
                    break
                end
            end
            
            # Check if we've closed a cycle in THIS trajectory
            if haskey(trajectory, x_new)
                cycle_start = trajectory[x_new]
                cycle_states = sort(collect(keys(trajectory)), by=s -> trajectory[s])
                cycle_states = cycle_states[(cycle_start + 1):end]
                if (cycle_states != [x_new])
                    push!(cycle_states, x_new)
                end
                
                # Check if this cycle matches an existing attractor
                attractor_id = find_equivalent_attractor(cycle_states, attractors)
                
                if attractor_id === nothing
                    # New attractor discovered
                    attractor_id = length(attractors) + 1
                    push!(attractors, cycle_states)
                    push!(attractor_counts, 1)
                    push!(attractor_times, [step])
                    push!(converged_flags, 1)
                else
                    # Existing attractor
                    attractor_counts[attractor_id] += 1
                    push!(attractor_times[attractor_id], step)
                end
                
                found_attractor = true
                break
            end
            
            trajectory[copy(x_new)] = step
            x = x_new
            
            if debug && step % 100 == 0
                println("Iteration $iter, Step $step")
            end
        end
        
        # Max steps reached without finding attractor
        if !found_attractor
            # Treat final state as an attractor (unconverged)
            attractor_id = find_equivalent_attractor([x], attractors)
            
            if attractor_id === nothing
                push!(attractors, [x])
                push!(attractor_counts, 1)
                push!(attractor_times, [max_steps])
                push!(converged_flags, 0)
            else
                attractor_counts[attractor_id] += 1
                push!(attractor_times[attractor_id], max_steps)
            end
        end
    end
    
    # Calculate statistics based on trajectory counts
    total_trajectories = length(initial_conditions)
    basin_proportions = attractor_counts ./ total_trajectories
    avg_times = [sum(times) / length(times) for times in attractor_times]
    
    
    df = DataFrame(
        states = attractors,
        basin_size = basin_proportions,
        time = avg_times,
        flag = converged_flags
    )
    
    return df
end

"""
    find_equivalent_attractor(cycle::Vector{Vector{Int}}, attractors::Vector{Vector{Vector{Int}}}) -> Union{Int, Nothing}

Check if a cycle matches an existing attractor (accounting for rotation).
"""
function find_equivalent_attractor(cycle::Vector{Vector{Int}}, attractors::Vector{Vector{Vector{Int}}})
    for (id, attractor) in enumerate(attractors)
        if cycles_equivalent(cycle, attractor)
            return id
        end
    end
    return nothing
end

"""
    cycles_equivalent(cycle1::Vector{Vector{Int}}, cycle2::Vector{Vector{Int}}) -> Bool

Check if two cycles are equivalent (same states, possibly rotated).
"""
function cycles_equivalent(cycle1::Vector{Vector{Int}}, cycle2::Vector{Vector{Int}})
    if length(cycle1) != length(cycle2)
        return false
    end
    
    n = length(cycle1)
    for offset in 0:(n-1)
        match = true
        for i in 1:n
            if cycle1[i] != cycle2[mod1(i + offset, n)]
                match = false
                break
            end
        end
        if match
            return true
        end
    end
    
    return false
end



# Example usage
"""
# Synchronous
F, I, N, degrees, variables, constants = getNodeFunctions("EMT_Switch.txt")

df_sync, state_dict = find_attractors_synchronous(
    F, N;
    nsim = 1000,
    exact = false,
    max_steps = 500,
    seed = 42
)

# Asynchronous
df_async = find_attractors_asynchronous(
    F, N;
    nsim = 10000,  # Need more samples for stochastic dynamics
    max_steps = 2000,  # Need more steps since only one node updates per step
    seed = 42
)
"""