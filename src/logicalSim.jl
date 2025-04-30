# boolean_network_sim.jl

"""
    update(F, I, N, X) -> Vector{Int}

Computes the updated state vector from the Boolean functions `F`,
input index list `I`, and current state vector `X`.

- `F`: List of Boolean functions (each a vector of length 2^k).
- `I`: List of regulator indices for each node.
- `N`: Number of nodes.
- `X`: Current state vector of size `N`.
"""
function update(F, I, N, X)
    Fx = zeros(Int, N)
    for i in 1:N
        Fx[i] = F[i][bin2dec(X[I[i]]) + 1]
    end
    return Fx
end

"""
    update_single_node(f, states_regulators) -> Int

Evaluates a single node's Boolean function `f` using the current
states of its regulators.
"""
function update_single_node(f, states_regulators)
    return f[bin2dec(states_regulators) + 1]
end

"""
    num_of_steady_states_asynchronous(F, I, N; kwargs...) -> Tuple

Runs an asynchronous simulation to find steady states in a Boolean network.

Arguments:
- `F`: Vector of Boolean functions (each a vector of Ints of length 2^k).
- `I`: Vector of regulator index vectors.
- `N`: Number of nodes.

Keyword arguments:
- `nsim`: Number of simulations to run (default: 500).
- `EXACT`: Whether to exhaustively simulate all 2^N states (default: false).
- `left_side_of_truth_table`: Used when `EXACT = true`.
- `initial_sample_points`: Optional specific starting points.
- `search_depth`: Max steps to simulate from each start (default: 50).
- `SEED`: RNG seed (default: -1, random).
- `DEBUG`: Whether to print debug info (default: true).

Returns a tuple:
    (steady_states, number_of_steady_states, basin_sizes,
     steady_state_dict, update_cache_dict, seed_used,
     sampled_points)
"""
function num_of_steady_states_asynchronous(
    F, I, N;
    nsim = 500,
    EXACT = false,
    left_side_of_truth_table = [],
    initial_sample_points = [],
    search_depth = 1000,
    SEED = -1,
    DEBUG = true,
    dictF = Dict{Tuple{Int, Int}, Int}()
)
    # Generate all 2^N binary input vectors if EXACT mode is enabled
    if EXACT && isempty(left_side_of_truth_table)
        left_side_of_truth_table = [reverse(collect(t)) for t in IterTools.product((0:1 for _ in 1:N)...)]
    end

    sampled_points = Int[]

    # Prevent using EXACT with initial sample points
    @assert isempty(initial_sample_points) || !EXACT "Sample points provided but with EXACT=true — ignored."

    # Set seed for reproducibility
    if SEED == -1
        SEED = rand(UInt32)
    end
    Random.seed!(SEED)

    # Dictionary to memoize function outputs for efficiency
    if length(dictF) == 0
        dictF = Dict{Tuple{Int, Int}, Int}()
    end

    # Data structures to store results
    steady_states = Int[]                  # list of final steady states (as decimal integers)
    basin_sizes = Int[]                    # basin size for each steady state
    steady_state_dict = Dict{Int, Int}()   # maps state -> (index in steady_states array - 1)
    timeVec = Int[]
    flagVec = Int[]
    frustrations = Float64[]

    iterations = EXACT ? 2^N : nsim        # how many initial conditions to try

    if !EXACT && isempty(initial_sample_points)
        # Generate random initial conditions
        initial_sample_points = [rand([0, 1], N) for i in 1:iterations]
    end
    counter = 0
    # Main simulation loop
    @inbounds for iteration in 1:iterations
        # Initialize a new state vector `x`
        if EXACT
            x = copy(left_side_of_truth_table[iteration])
            xbin = iteration - 1
        else
            x = copy(initial_sample_points[iteration])
            xbin = bin2dec(x)
        end

        if DEBUG
            println(iteration," ", -1, -1, false, xbin, x)
        end
        HAS_ATTRACTOR = Ref(false)
        # Run asynchronous updates for at most `search_depth` steps
        for jj in 1:search_depth
            FOUND_NEW_STATE = false
            # This is true if all possible updates of the state are not known


            # If current state already known as steady state, assign basin and break
            index_ss = get(steady_state_dict, xbin, nothing)
            if !isnothing(index_ss)
                # Found a steady state
                basin_sizes[index_ss + 1] += 1
                break
            else
                # Try one asynchronous update round in a random order
                update_order_to_try = randperm(N)
                for i in update_order_to_try # Until a steady state or a new state is found
                    # Try to retrieve cached result of updating node i from this state
                    fxbin = get(dictF, (xbin, i), nothing)
                    # If fxbin is nothing, we don't know the update, hence update. 
                    # check if xbin is unvarying at i. if not, the state gets updated and we found a new state.
                    if fxbin !== nothing
                        if fxbin != xbin
                            # update is known, but not the same as xbin. So x[i] flips on update
                            FOUND_NEW_STATE = true
                            x[i] = 1 - x[i]
                        end
                    else
                        # Apply node update rule
                        fx_i = update_single_node(F[i], x[I[i]])
                        if fx_i == x[i]
                            fxbin = xbin
                        else
                            # if 0 gets updated to 1, add 2^(N-1) else subtract
                            fxbin = xbin + (fx_i - x[i])*2^(N - i)
                            FOUND_NEW_STATE = true
                        end
                        x[i] = fx_i
                        # Cache result
                        dictF[(xbin, i)] = fxbin
                    end
                    # x has been updated
                    if FOUND_NEW_STATE
                        xbin = fxbin
                        break
                    end
                    # xbin has been updated if a new state was found. 
                    #That is, there is an i for which fx_i != x[i], and therefore xbin was not a steady state
                end

                if DEBUG
                    println(iteration, jj, i, FOUND_NEW_STATE, xbin, x)
                end
            end
            # if at any jj, we find a steady state, FOUND_NEW_STATE remains false, therefore we get out of the loop and record the jj as time
            if !FOUND_NEW_STATE
                HAS_ATTRACTOR[] = true
                
                # No further update — reached a steady state
                index_ss = get(steady_state_dict, xbin, nothing)
                if index_ss !== nothing
                    basin_sizes[index_ss + 1] += 1
                    timeVec[index_ss + 1] = timeVec[index_ss + 1] + jj
                else
                    steady_state_dict[xbin] = length(steady_states)
                    push!(steady_states, xbin)
                    push!(basin_sizes, 1)
                    push!(flagVec, 1)
                    push!(timeVec, jj)
                end
                continue
            end
        end
        if !HAS_ATTRACTOR[]
            index_ss = get(steady_state_dict, xbin, nothing)
            if index_ss !== nothing
                basin_sizes[index_ss + 1] += 1
                timeVec[index_ss + 1] = timeVec[index_ss + 1] + search_depth
            else
                steady_state_dict[xbin] = length(steady_states)
                push!(steady_states, xbin)
                push!(basin_sizes, 1)
                push!(flagVec, 0)
                push!(timeVec, search_depth)
            end
        end
        if DEBUG
            println()
        end
    end
    basin_sizes = basin_sizes./sum(basin_sizes)
    ssDf = DataFrame(
        steady_state = steady_states,
        basin_size = basin_sizes,
        time_to_reach = timeVec./sum(basin_sizes),
        flag = flagVec
    )
    # Final check: Warn if some runs didn’t converge to a state
    return ssDf, dictF
end

