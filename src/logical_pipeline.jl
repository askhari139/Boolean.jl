# complete_boolean_pipeline.jl
# End-to-end pipeline: Rules → Truth Tables → Attractors → Perturbations

using DataFrames
using Random
using CSV

"""
    f_from_expression(expr::String) -> Tuple{Vector{Int}, Vector{String}}

Parse a Boolean expression and return (truth_table, variable_names).
Uses the tokenize/postfix/DNF pipeline.
"""
function f_from_expression(expr::String)
    tokens = tokenize_expression(expr)
    postfix = infix_to_postfix(tokens)
    dnf = postfix_to_dnf(postfix)
    tt_df = dnf_to_truthtable(dnf)
    
    # Extract variable names (all columns except "Output")
    vars = [name for name in names(tt_df) if name != "Output"]
    
    # Extract truth table values
    truth_table = Int.(tt_df.Output)
    
    return truth_table, vars
end

"""
    parse_boolean_network(rules_file::String) -> Tuple

Complete parsing pipeline: reads Boolean rules and returns simulation-ready format.

Arguments:
- `rules_file`: Path to file with Boolean rules (format: "NodeName = Boolean Expression")

Returns:
- `F`: Vector of truth tables for each node
- `I`: Vector of regulator index vectors
- `N`: Number of nodes
- `node_names`: Vector of node names
- `node_to_idx`: Dictionary mapping node names to indices
"""
function parse_boolean_network(rules_file::String)
    # Use existing parser
    truth_table_dict = f_from_file(rules_file)
    
    # Create node name to index mapping
    node_names = sort(collect(keys(truth_table_dict)))
    N = length(node_names)
    node_to_idx = Dict(name => i for (i, name) in enumerate(node_names))
    
    # Initialize F and I
    F = Vector{Vector{Int}}(undef, N)
    I = Vector{Vector{Int}}(undef, N)
    
    # Fill F and I for each node
    for (i, node_name) in enumerate(node_names)
        truth_table, regulators = truth_table_dict[node_name]
        
        # Store truth table
        F[i] = truth_table
        
        # Convert regulator names to indices
        if isempty(regulators)
            # Constitutive node - self-regulation convention
            I[i] = [i]
        else
            regulator_indices = [node_to_idx[reg] for reg in regulators]
            I[i] = regulator_indices
        end
    end
    
    return F, I, N, node_names, node_to_idx
end

"""
    find_attractors_synchronous(F, I, N; kwargs...) -> DataFrame

Find attractors using synchronous (parallel) updates.
All nodes update simultaneously at each time step.
"""
function find_attractors_synchronous(
    F, I, N;
    nsim = 500,
    EXACT = false,
    search_depth = 1000,
    SEED = -1,
    DEBUG = false
)
    if SEED == -1
        SEED = rand(UInt32)
    end
    Random.seed!(SEED)
    
    # Data structures
    attractors = Dict{Vector{Int}, Int}()  # state_vector => attractor_id
    attractor_states = Vector{Vector{Int}}()  # List of attractor state sequences
    basin_sizes = Int[]
    
    # Generate initial conditions
    if EXACT
        left_side = [reverse(collect(t)) for t in Iterators.product((0:1 for _ in 1:N)...)]
        iterations = 2^N
        initial_states = left_side
    else
        iterations = nsim
        initial_states = [rand([0, 1], N) for _ in 1:iterations]
    end
    
    # Simulate from each initial condition
    for iteration in 1:iterations
        x = copy(initial_states[iteration])
        trajectory = [copy(x)]
        
        for step in 1:search_depth
            # Synchronous update: compute all next states simultaneously
            x_next = zeros(Int, N)
            for i in 1:N
                x_next[i] = update_single_node(F[i], x[I[i]])
            end
            
            x = x_next
            push!(trajectory, copy(x))
            
            # Check if we've seen this state in trajectory (cycle detection)
            cycle_start = findfirst(==(x), trajectory[1:end-1])
            if !isnothing(cycle_start)
                # Found a cycle
                cycle = trajectory[cycle_start:end-1]
                
                # Check if this attractor already exists
                attractor_id = nothing
                for (existing_cycle, id) in attractors
                    if length(existing_cycle) == length(cycle) && 
                       any(circshift(cycle, k) == existing_cycle for k in 0:length(cycle)-1)
                        attractor_id = id
                        break
                    end
                end
                
                if isnothing(attractor_id)
                    # New attractor
                    attractor_id = length(attractor_states) + 1
                    push!(attractor_states, cycle)
                    push!(basin_sizes, 1)
                    for state in cycle
                        attractors[copy(state)] = attractor_id
                    end
                else
                    basin_sizes[attractor_id] += 1
                end
                
                break
            end
        end
    end
    
    # Build results DataFrame
    results = DataFrame(
        attractor_id = Int[],
        cycle_length = Int[],
        representative_state = Int[],  # First state as decimal
        basin_size = Float64[],
        is_fixed_point = Bool[]
    )
    
    for (i, cycle) in enumerate(attractor_states)
        rep_state = bin2dec(cycle[1])
        push!(results, (
            attractor_id = i,
            cycle_length = length(cycle),
            representative_state = rep_state,
            basin_size = basin_sizes[i] / sum(basin_sizes),
            is_fixed_point = length(cycle) == 1
        ))
    end
    
    return results, attractor_states
end

"""
    analyze_network_complete(rules_file::String; kwargs...) -> Dict

Complete analysis pipeline from rules file to attractors.

Returns a dictionary with:
- `network_info`: Basic network information
- `asynchronous_attractors`: DataFrame of async attractors
- `synchronous_attractors`: DataFrame of sync attractors
- `F`, `I`, `N`, `node_names`: Network specification
"""
function analyze_network_complete(
    rules_file::String;
    nsim_async = 500,
    nsim_sync = 500,
    EXACT = false,
    search_depth = 1000,
    SEED = -1,
    DEBUG = false
)
    println("="^60)
    println("Complete Boolean Network Analysis Pipeline")
    println("="^60)
    
    # Step 1: Parse network
    println("\n[1/3] Parsing Boolean rules from file...")
    F, I, N, node_names, node_to_idx = parse_boolean_network(rules_file)
    println("  ✓ Parsed $N nodes")
    
    # Network statistics
    avg_degree = mean(length(regs) for regs in I)
    max_degree = maximum(length(regs) for regs in I)
    println("  ✓ Average in-degree: $(round(avg_degree, digits=2))")
    println("  ✓ Maximum in-degree: $max_degree")
    
    # Step 2: Find asynchronous attractors
    println("\n[2/3] Finding asynchronous attractors...")
    async_df, async_cache = num_of_steady_states_asynchronous(
        F, I, N,
        nsim = nsim_async,
        EXACT = EXACT,
        search_depth = search_depth,
        SEED = SEED,
        DEBUG = DEBUG
    )
    println("  ✓ Found $(nrow(async_df)) asynchronous attractor(s)")
    
    # Step 3: Find synchronous attractors
    println("\n[3/3] Finding synchronous attractors...")
    sync_df, sync_attractors = find_attractors_synchronous(
        F, I, N,
        nsim = nsim_sync,
        EXACT = EXACT,
        search_depth = search_depth,
        SEED = SEED,
        DEBUG = DEBUG
    )
    println("  ✓ Found $(nrow(sync_df)) synchronous attractor(s)")
    
    # Summary
    println("\n" * "="^60)
    println("Summary")
    println("="^60)
    println("Asynchronous steady states: $(nrow(async_df))")
    for row in eachrow(async_df)
        state_vec = dec2bin(row.steady_state, N)
        println("  State $(row.steady_state): $state_vec")
        println("    Basin: $(round(row.basin_size * 100, digits=1))%")
        println("    Converged: $(row.flag == 1 ? "Yes" : "No")")
    end
    
    println("\nSynchronous attractors: $(nrow(sync_df))")
    for row in eachrow(sync_df)
        println("  Attractor $(row.attractor_id):")
        println("    Cycle length: $(row.cycle_length)")
        println("    Fixed point: $(row.is_fixed_point)")
        println("    Basin: $(round(row.basin_size * 100, digits=1))%")
    end
    
    return Dict(
        "network_info" => Dict(
            "N" => N,
            "node_names" => node_names,
            "avg_degree" => avg_degree,
            "max_degree" => max_degree
        ),
        "asynchronous_attractors" => async_df,
        "synchronous_attractors" => sync_df,
        "synchronous_attractor_cycles" => sync_attractors,
        "F" => F,
        "I" => I,
        "N" => N,
        "node_names" => node_names,
        "node_to_idx" => node_to_idx,
        "cache" => async_cache
    )
end

"""
    perturb_attractors_complete(
        results::Dict,
        perturbation_sizes = [1, 2];
        kwargs...
    ) -> DataFrame

Complete perturbation analysis on all attractors found in the network.

Arguments:
- `results`: Output from analyze_network_complete()
- `perturbation_sizes`: Vector of perturbation sizes to test
- `ntrials_per_size`: Number of trials per perturbation size
- `max_steps`: Maximum simulation steps
- `update_mode`: :asynchronous or :synchronous

Returns DataFrame with perturbation results for all attractors.
"""
function perturb_attractors_complete(
    results::Dict;
    perturbation_sizes = [1, 2],
    ntrials_per_size = 100,
    max_steps = 1000,
    update_mode = :asynchronous,
    SEED = -1
)
    println("\n" * "="^60)
    println("Perturbation Analysis")
    println("="^60)
    
    F = results["F"]
    I = results["I"]
    N = results["N"]
    node_names = results["node_names"]
    
    if SEED == -1
        SEED = rand(UInt32)
    end
    
    # Determine which attractors to perturb based on update mode
    if update_mode == :asynchronous
        attractors_df = results["asynchronous_attractors"]
        attractor_states = [row.steady_state for row in eachrow(attractors_df)]
        println("Analyzing $(length(attractor_states)) asynchronous steady states")
    else
        attractors_df = results["synchronous_attractors"]
        attractor_cycles = results["synchronous_attractor_cycles"]
        # For synchronous, perturb the first state of each cycle
        attractor_states = [bin2dec(cycle[1]) for cycle in attractor_cycles]
        println("Analyzing $(length(attractor_states)) synchronous attractors")
    end
    
    # Collect all perturbation results
    all_results = DataFrame(
        original_attractor = Int[],
        perturbation_size = Int[],
        nodes_flipped = Vector{Int}[],
        nodes_flipped_names = Vector{String}[],
        final_attractor = Int[],
        steps_to_steady = Int[],
        converged = Bool[],
        returned_to_original = Bool[]
    )
    
    # Analyze each attractor
    for (idx, attractor_decimal) in enumerate(attractor_states)
        println("\n--- Attractor $idx (State: $attractor_decimal) ---")
        
        pert_results = analyze_perturbation_response(
            F, I, N, attractor_decimal,
            perturbation_sizes = perturbation_sizes,
            ntrials_per_size = ntrials_per_size,
            max_steps = max_steps,
            SEED = SEED + idx  # Different seed for each attractor
        )
        
        # Add attractor ID and node names
        for row in eachrow(pert_results)
            node_names_flipped = [node_names[i] for i in row.nodes_flipped]
            
            push!(all_results, (
                original_attractor = attractor_decimal,
                perturbation_size = row.perturbation_size,
                nodes_flipped = row.nodes_flipped,
                nodes_flipped_names = node_names_flipped,
                final_attractor = row.final_state,
                steps_to_steady = row.steps_to_steady,
                converged = row.converged,
                returned_to_original = row.returned_to_original
            ))
        end
        
        # Summary for this attractor
        return_rate = mean(pert_results.returned_to_original)
        avg_steps = mean(pert_results.steps_to_steady)
        
        println("  Return probability: $(round(return_rate * 100, digits=1))%")
        println("  Average recovery time: $(round(avg_steps, digits=1)) steps")
        
        # Show transition destinations
        transitions = combine(groupby(pert_results, :final_state), 
                              nrow => :count)
        sort!(transitions, :count, rev=true)
        
        println("  Transition destinations:")
        for trans_row in eachrow(transitions)
            prob = trans_row.count / nrow(pert_results) * 100
            dest_state = dec2bin(trans_row.final_state, N)
            println("    → State $(trans_row.final_state) ($dest_state): $(round(prob, digits=1))%")
        end
    end
    
    # Overall stability comparison
    println("\n" * "="^60)
    println("Stability Comparison Across Attractors")
    println("="^60)
    
    stability_summary = combine(
        groupby(all_results, [:original_attractor, :perturbation_size]),
        :returned_to_original => mean => :return_probability,
        :steps_to_steady => mean => :avg_recovery_time,
        nrow => :n_trials
    )
    
    println(stability_summary)
    
    # Find most and least stable
    for pert_size in perturbation_sizes
        subset = stability_summary[stability_summary.perturbation_size .== pert_size, :]
        if nrow(subset) > 0
            most_stable_idx = argmax(subset.return_probability)
            least_stable_idx = argmin(subset.return_probability)
            
            println("\nFor $(pert_size)-node perturbations:")
            println("  Most stable: Attractor $(subset.original_attractor[most_stable_idx])")
            println("    Return probability: $(round(subset.return_probability[most_stable_idx] * 100, digits=1))%")
            println("  Least stable: Attractor $(subset.original_attractor[least_stable_idx])")
            println("    Return probability: $(round(subset.return_probability[least_stable_idx] * 100, digits=1))%")
        end
    end
    
    return all_results, stability_summary
end

"""
    full_network_analysis(
        rules_file::String;
        kwargs...
    ) -> Dict

Complete end-to-end analysis: parse → find attractors → perturb → analyze.

This is the main function to call for a complete analysis.
"""
function full_network_analysis(
    rules_file::String;
    nsim = 500,
    EXACT = false,
    search_depth = 1000,
    perturbation_sizes = [1, 2],
    ntrials_per_perturbation = 100,
    analyze_both_modes = false,
    output_dir = nothing,
    SEED = -1
)
    println("╔" * "="^58 * "╗")
    println("║" * " "^15 * "BOOLEAN NETWORK ANALYSIS" * " "^19 * "║")
    println("╚" * "="^58 * "╝")
    println("\nInput file: $rules_file")
    println("Simulations: $nsim")
    println("EXACT mode: $EXACT")
    println()
    
    # Step 1: Analyze network and find attractors
    results = analyze_network_complete(
        rules_file,
        nsim_async = nsim,
        nsim_sync = nsim,
        EXACT = EXACT,
        search_depth = search_depth,
        SEED = SEED,
        DEBUG = false
    )
    
    # Step 2: Perturbation analysis on asynchronous attractors
    println("\n" * "╔" * "="^58 * "╗")
    println("║" * " "^12 * "ASYNCHRONOUS PERTURBATIONS" * " "^20 * "║")
    println("╚" * "="^58 * "╝")
    
    async_pert_results, async_stability = perturb_attractors_complete(
        results,
        perturbation_sizes = perturbation_sizes,
        ntrials_per_size = ntrials_per_perturbation,
        update_mode = :asynchronous,
        SEED = SEED
    )
    
    # Step 3: Optionally analyze synchronous mode
    sync_pert_results = nothing
    sync_stability = nothing
    
    if analyze_both_modes
        println("\n" * "╔" * "="^58 * "╗")
        println("║" * " "^13 * "SYNCHRONOUS PERTURBATIONS" * " "^20 * "║")
        println("╚" * "="^58 * "╝")
        
        sync_pert_results, sync_stability = perturb_attractors_complete(
            results,
            perturbation_sizes = perturbation_sizes,
            ntrials_per_size = ntrials_per_perturbation,
            update_mode = :synchronous,
            SEED = SEED
        )
    end
    
    # Step 4: Save results if output directory specified
    if !isnothing(output_dir)
        mkpath(output_dir)
        
        base_name = splitext(basename(rules_file))[1]
        
        # Save attractor information
        CSV.write(
            joinpath(output_dir, "$(base_name)_async_attractors.csv"),
            results["asynchronous_attractors"]
        )
        
        CSV.write(
            joinpath(output_dir, "$(base_name)_sync_attractors.csv"),
            results["synchronous_attractors"]
        )
        
        # Save perturbation results
        CSV.write(
            joinpath(output_dir, "$(base_name)_async_perturbations.csv"),
            async_pert_results
        )
        
        CSV.write(
            joinpath(output_dir, "$(base_name)_async_stability.csv"),
            async_stability
        )
        
        if !isnothing(sync_pert_results)
            CSV.write(
                joinpath(output_dir, "$(base_name)_sync_perturbations.csv"),
                sync_pert_results
            )
            CSV.write(
                joinpath(output_dir, "$(base_name)_sync_stability.csv"),
                sync_stability
            )
        end
        
        println("\n✓ Results saved to: $output_dir")
    end
    
    # Return everything
    return Dict(
        "network" => results,
        "async_perturbations" => async_pert_results,
        "async_stability" => async_stability,
        "sync_perturbations" => sync_pert_results,
        "sync_stability" => sync_stability
    )
end

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

"""
Example usage:

# Simple analysis
results = full_network_analysis(
    "my_rules.txt",
    nsim = 500,
    EXACT = false,
    perturbation_sizes = [1, 2],
    ntrials_per_perturbation = 100,
    output_dir = "results"
)

# Access specific results
async_attractors = results["network"]["asynchronous_attractors"]
perturbation_data = results["async_perturbations"]
stability_metrics = results["async_stability"]

# For small networks, use EXACT mode
results_exact = full_network_analysis(
    "small_network.txt",
    EXACT = true,
    perturbation_sizes = [1],
    output_dir = "results_exact"
)
"""
