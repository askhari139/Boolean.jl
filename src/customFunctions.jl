### Edge weights

function edgeWeightPert(topoFile::String; nPerts::Int=10000, nInit::Int64=10000, nIter::Int64=1000, 
    mode::String="Async", stateRep::Int64=-1, reps::Int = 3, csv::Bool=false, 
    types::Array{Int, 1} = [0],
    minWeight::Float64=0.0, maxWeight::Float64=1.0,
    newFolder::Bool = true)
    updMat, nodes = topo2interaction(topoFile)
    nzId = enumerate(findall(updMat.!=0))
    edges = [join([nodes[j[1]], "_", nodes[j[2]]]) for (i,j) in nzId]
    nZ = length(nzId)
    nRand = nZ*nPerts
    rands = rand(nRand)
    rands = minWeight .+ rands*(maxWeight-minWeight)
    rands = reshape(rands, nPerts, nZ)
    randFold = replace(topoFile, ".topo" => "_rand")
    d1 = pwd()
    mkpath(randFold)
    cpPath = join([randFold, "/", topoFile])
    cp(topoFile, cpPath, force = true)
    cd(randFold)
    if newFolder
        mkpath(join([minWeight, maxWeight], "_"))
        cpPath = join([minWeight, "_", maxWeight, "/", topoFile])
        cp(topoFile, cpPath, force = true)
        cd(join([minWeight, "_", maxWeight]))
        colNames = [Symbol("rand", i) for i in 1:nPerts]
        append!(colNames, [:edge])
        df = DataFrame(hcat(rands', edges), colNames)
        CSV.write(replace(topoFile, ".topo" => "_edgeWeightPert.dat"), df; delim = " ")

    end
    nodesName = join([replace(topoFile, ".topo" => ""), "_nodes.txt"])
    update_matrix,Nodes = topo2interaction(topoFile)
    io = open(nodesName, "w")
    for i in Nodes
        println(io, i)
    end
    close(io);
    p = Progress(nPerts)
    Threads.@threads for i in 1:nPerts
        # println(string(i))
        bmodel_reps(topoFile; nInit = nInit, nIter = nIter, mode = mode,
        stateRep = stateRep, randSim=true, root = string(i), 
        randVec = rands[i,:], types = types, reps = reps)
    end
    finish!(p)
    
    cd(d1)

end

function weightedTopoSim(topoFiles::Vector{String}; nInit::Int64=10000, 
    nIter::Int64=1000, mode::String="Async", stateRep::Int64=-1, reps::Int = 3)
    p = Progress(nPerts)
    Threads.@threads for topoFile in topoFiles
        # println(string(i))
        bmodel_reps(topoFile; nInit = nInit, nIter = nIter, mode = mode,
        stateRep = stateRep, randSim=true, root = string(i), 
        randVec = rands[i,:], types = types, reps = reps)
    end
    finish!(p)
    cd(d1)
end



function contWeightPert(topoFile::String; nInit::Int64=1000,
    nIter::Int64=100000, mode::String="Async", stateRep::Int64=-1, 
    noise::Float64=0.01, steadyStates::Bool=true, topN::Int=10,
    frequency::Int=1)
    
    # Load network topology
    updMat, nodes = topo2interaction(topoFile)
    
    # Run simulation with proper weight function
    weightFunc = defaultWeightsFunction(noise)
    y1 = @elapsed stateMat, states = asyncRandCont(updMat, nInit, nIter, stateRep; 
        noise=noise, steadyStates=steadyStates, topN=topN, frequency=frequency)
    println("Async rand took $y1")
    nStates = length(states)
    
    # Calculate Mean Residence Time (MRT) for each state - VECTORIZED
    stateFreqs = zeros(Float64, nInit, nStates)
    # stateFreqs = Dict{Int, Vector{Float64}}()
    Threads.@threads for i in 1:nInit
        counts = Dict{Int, Int}()
        for j in 1:nIter
            stateID = stateMat[i,j]
            counts[stateID] = get(counts, stateID, 0) + 1
        end
        for (stateID, count) in counts
            stateFreqs[i, stateID] = count / nIter
        end
    end
    
    MRT = vec(mean(stateFreqs, dims=1)) # column-wise means - each state across all trajectories
    MRTsd = vec(std(stateFreqs, dims=1)) # column-wise std
    
    # Count switching events per trajectory - VECTORIZED
    switchingEvents = zeros(Int, nInit)
    for i in 1:nInit
        switchingEvents[i] = sum(stateMat[i, 2:end] .!= stateMat[i, 1:end-1])
    end
    
    # Save results
    baseName = replace(topoFile, ".topo" => "")
    
    # Save trajectory matrix
    trajectoryDF = DataFrame(stateMat, :auto)
    rename!(trajectoryDF, ["t$i" for i in 1:nIter])
    CSV.write("$(baseName)_contWeightPert.csv", trajectoryDF)
    
    # Save state statistics
    stateStatsDF = DataFrame(
        ID = 1:nStates,
        State = states,
        MRT = MRT,
        MRTsd = MRTsd
    )
    sort!(stateStatsDF, :MRT, rev=true)
    CSV.write("$(baseName)_contWeightPert_states.csv", stateStatsDF)
    
    # Extract initial states from first column of stateMat
    initialStateIDs = stateMat[:, 1]
    initialStates = states[initialStateIDs]
    
    # Save trajectory-level statistics
    trajectoryStatsDF = DataFrame(
        TrajectoryID = 1:nInit,
        InitialStateID = initialStateIDs,
        InitialState = initialStates,
        SwitchingEvents = switchingEvents
    )
    CSV.write("$(baseName)_contWeightPert_trajectories.csv", trajectoryStatsDF)
    
    # If using steady states, add per-initial-condition statistics
    if steadyStates
        # Group by initial state and calculate statistics
        initialConditionStats = combine(groupby(trajectoryStatsDF, [:InitialStateID, :InitialState])) do df
            trajIDs = df.TrajectoryID
            
            # Calculate MRT for each state, averaged over trajectories from this initial condition
            stateMRTs = [mean(stateFreqs[trajIDs, stateID]) for stateID in 1:nStates]
            
            # Find top 5 most visited states
            topIndices = sortperm(stateMRTs, rev=true)[1:min(5, nStates)]
            
            DataFrame(
                nTrajectories = nrow(df),
                MeanSwitchingEvents = mean(df.SwitchingEvents),
                StdSwitchingEvents = std(df.SwitchingEvents),
                MedianSwitchingEvents = median(df.SwitchingEvents),
                TopStates = join(states[topIndices], "; "),
                TopStatesMRT = join(round.(stateMRTs[topIndices], digits=4), "; ")
            )
        end
        
        CSV.write("$(baseName)_contWeightPert_byInitialState.csv", initialConditionStats)
    end
    
    return stateStatsDF, trajectoryStatsDF
end

function getSSListRand(topoFile::String;
    minWeight::Float64=0.0, maxWeight::Float64=1.0, nPerts::Int=10000)
    updMat, nodes = topo2interaction(topoFile)
    nzId = enumerate(findall(updMat.!=0))
    edges = [join([nodes[j[1]], "_", nodes[j[2]]]) for (i,j) in nzId]
    nZ = length(nzId)
    nRand = nZ*nPerts
    rands = rand(nRand)
    rands = minWeight .+ rands*(maxWeight-minWeight)
    rands = reshape(rands, nPerts, nZ)
    n_nodes = size(updMat,1)
    nStates = min(2^n_nodes, 100000)
    stateVec = [-1.0,1.0]
    if (n_nodes < 20)
        states = listStates(n_nodes, stateVec)
    else
        states = [rand(stateVec, n_nodes) for i in 1:nStates]
    end
    # print(states)
    finVecAll = []
    @showprogress for i in 1:nPerts
        randVec = rands[i,:]
        update_matrix = copy(updMat)
        nzId = enumerate(findall(update_matrix.!=0))
        update_matrix = [Float64(i) for i in update_matrix]
        for (i,j) in nzId
            update_matrix[j] = update_matrix[j]*randVec[i]
        end
        # print(update_matrix)
        # return
        finVec = []
        for state in states
            s1 = float(sign.(update_matrix*state))
            s1 = [s1[i] == 0 ? state[i] : s1[i] for i in 1:n_nodes]
            if s1 == state
                fin = join(Int.(replace(x -> x == -1 ? 0 : 1, state)))
                push!(finVec, fin)
            end
        end
        finVec = join(sort(finVec), " ")
        push!(finVecAll, finVec)
    end
    ssFile = replace(topoFile, ".topo" => "_ssRand")
    ssFile = join([ssFile, "_", minWeight, "_", maxWeight, ".dat"])
    io = open(ssFile, "w")
    for i in finVecAll
        println(io, i)
    end
    close(io);
end

### state transition graph - incomplete
function stg(topoFile::String, mode::String="Async")
    print(topoFile)
    update_matrix,Nodes = topo2interaction(topoFile)
    if mode == "Async"
        stg = async_stg(update_matrix)
    else
        stg = sync_stg(update_matrix)
    end
    CSV.write(join([replace(topoFile, ".topo" => ""), "_stg.csv"]), stg)
    return stg,Nodes
end

### single Node turnOff 
function scanNodeTurnOff(topoFile::String; nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, reps::Int = 3, init::Bool=false, 
    root::String="", nLevels = 2, progressMeter::Bool = false, 
    getData::Bool = false, tSet::Array{Int, 1} = Int64[], getPrevious::Bool = false)
    update_matrix,Nodes = topo2interaction(topoFile)
    n_nodes = size(update_matrix,1)
    if getPrevious
        if length(tSet) == 0
        ### Normal simulation
            dfs = bmodel_reps(topoFile; nInit = nInit, nIter = nIter, 
                    mode = mode, stateRep = stateRep, reps = reps, init = init, nLevels = nLevels,
                    root = root, shubham = true, vaibhav = false, write = false, getData = true)
            if init
                init, finFlag = dfs
                # add a turnOffNode column to init and finFlag with the value "None"
                init[!, "turnOffNode"] .= "None"
                initList = [init]
            else
                finFlag = dfs
            end
            finFlag[!, "turnOffNode"] .= "None"
            finFlagList = [finFlag]
        else
            dfs = bmodel_reps(topoFile; nInit = nInit, nIter = nIter, 
                mode = mode, stateRep = stateRep, reps = reps, init = init, nLevels = nLevels,
                root = root, shubham = true, vaibhav = true, write = false, getData = true, 
                turnOffNodes = tSet)
            if init
                init, finFlag = dfs
                initList = [init]
            else
                finFlag = dfs
            end
            finFlagList = [finFlag]
        end
    else
        initList = []
        finFlagList = []
    end
    
    Threads.@threads for i in 1:n_nodes
        if (i == n_nodes + 1)
            turnOffNodes = Int[]
        elseif i in tSet
            continue
        else
            turnOffNodes = vcat(tSet, i)
        end
        dfs = bmodel_reps(topoFile; nInit = nInit, nIter = nIter, 
            mode = mode, stateRep = stateRep, reps = reps, init = init, nLevels = nLevels,
            root = root, shubham = true, vaibhav = true, 
            turnOffNodes = turnOffNodes, write = false, getData = true)
        if init
            finFlag, init = dfs
            push!(initList, init)
        else
            finFlag = dfs
        end
        push!(finFlagList, finFlag)
    end
        if init
            initDf = vcat(initList...)
        end
        finFlagDf = vcat(finFlagList...)
        if init
            CSV.write(replace(topoFile, ".topo" => "_scanNode_turnOff_initFinFlagFreq.csv"), initDf)
        end
        CSV.write(replace(topoFile, ".topo" => "_scanNode_turnOff_finFlagFreq.csv"), finFlagDf)
    if getData
        if init
            return initDf, finFlagDf
        else
            return finFlagDf
        end
    end
end

"""
    format_attractor_states(states::Vector{Vector{Int}}) -> String

Format attractor states as a readable string.
"""
function format_attractor_states(states::Vector{Vector{Int}})
    state_strings = [join(state, "_") for state in states]
    return join(state_strings, ";")
end

"""
    simulate_network_folder(folder::String; 
                            reps::Int=10, 
                            output_csv::String="results.csv", 
                            parallel::Bool=false,
                            # text_to_BN parameters
                            separator_var_func::String="=", 
                            original_not::String="NOT", 
                            original_and::String="AND", 
                            original_or::String="OR", 
                            new_not::String=" not ", 
                            new_and::String=" and ", 
                            new_or::String=" or ", 
                            max_degree::Int=15, 
                            TREATMENT_OF_CONSTANTS::Int=1, 
                            max_N::Int=10000,
                            # num_of_steady_states_asynchronous parameters
                            max_iterations::Int=200, 
                            nstarts::Int=100, 
                            cutoff::Int=10, 
                            report_freq::Int=0, 
                            sample_frac::Float64=1.0, 
                            seed::Union{Int,Nothing}=nothing)

Simulates all Boolean networks in the specified folder. Each network is simulated `reps` times.

# Arguments
- `folder`: Folder path containing text-based Boolean network definitions.
- `reps`: Number of repetitions per network (default 10).
- `output_csv`: Output CSV file path (default "results.csv").
- `parallel`: If true, enables multithreaded execution across networks.
- `text_to_BN` arguments: Customize Boolean logic parsing and constants treatment.
- `num_of_steady_states_asynchronous` arguments: Customize convergence behavior and initial condition sampling.

# Output CSV columns
- `File`: Network filename.
- `State`: Unique steady state (as a string).
- `Frequency`: Average basin size (fraction of ICs converging to that state).
- `Time`: Average convergence time.
- `Flag`: 1 if valid steady state.
- `Frustration`: Network frustration score.

"""
function simulate_network_logical(
    rules_file::String;
    folder::String = ".",
    n_replicates::Int = 3,
    update_mode::String = "synchronous",
    max_steps::Int = 1000,
    n_initial_conditions::Int = 100000,
    exact::Bool = false,
    seed::Int = -1
)
    @assert update_mode in ["synchronous", "asynchronous"] "Invalid update_mode"
    
    # Load network
    rules_path = joinpath(folder, rules_file)
    F, I, N, degrees, variables, constants = getNodeFunctions(rules_path)
    
    # Adjust parameters based on update mode
    if update_mode == "asynchronous"
        exact = false
        n_initial_conditions = max(2^N, n_initial_conditions)
        max_steps = max(max_steps, N * 100)  # Need more steps for async
    else
        if exact && N > 17
            @warn "Network too large for exact mode (N=$N). Switching to sampling."
            exact = false
        end
    end
    
    # Load topology for frustration
    interaction_indices = nothing
    node_names = nothing
    
    # if calculate_frustration
    topo_file = replace(rules_file, ".txt" => ".topo")
    topo_path = joinpath(folder, topo_file)
    
    if !isfile(topo_path)
        boolean_to_topo(rules_path)
    end
    update_matrix, node_names = topo2interaction(topo_path)
    interaction_indices = findnz(sparse(update_matrix))
    # end
    
    # Run replicates
    println("Running $n_replicates replicates with $update_mode update mode...")
    all_results = DataFrame[]
    
    for rep in 1:n_replicates
        println("  Replicate $rep/$n_replicates")
        
        if update_mode == "synchronous"
            result_df, _ = Boolean.find_attractors_synchronous(
                F, N;
                nsim = n_initial_conditions,
                exact = exact,
                max_steps = max_steps,
                seed = seed,
                debug = false
            )
        else  # asynchronous
            result_df = find_attractors_asynchronous(
                F, N;
                nsim = n_initial_conditions,
                max_steps = max_steps,
                seed = seed,
                debug = false
            )
        end
        
        push!(all_results, result_df)
    end
    
    # Combine results
    combined_df = combine_replicate_results(all_results)
    
    # Add frustration
    combined_df = add_frustration_scores(
        combined_df, N, variables, constants, 
        node_names, interaction_indices
    )
    combined_df.states = format_attractor_states.(combined_df.states)
    
    sort!(combined_df, :basin_size_mean, rev=true)
    
    # Save results
    output_file = replace(rules_file, ".txt" => "_attractors.csv")
    output_path = joinpath(folder, output_file)
    CSV.write(output_path, combined_df)
    
    nodesName = replace(rules_file, ".txt" => "_nodes.txt")
    io = open(nodesName, "w")
    for i in node_names
        println(io, i)
    end
    close(io);

    println("Results saved to: $output_path")
    
    return combined_df
end
"""
    combine_replicate_results(results::Vector{DataFrame}) -> DataFrame

Combine results from multiple simulation replicates, computing mean and SD.
"""
function combine_replicate_results(results::Vector{DataFrame})
    # Join all dataframes by attractor states and flag
    combined = results[1]
    
    for i in 2:length(results)
        combined = outerjoin(
            combined, 
            results[i], 
            on = [:states, :flag], 
            makeunique = true
        )
    end
    
    # Calculate mean and SD for basin sizes
    basin_cols = [col for col in names(combined) if startswith(string(col), "basin_size")]
    
    combined.basin_size_mean = [
        mean(skipmissing([row[col] for col in basin_cols]))
        for row in eachrow(combined)
    ]
    
    combined.basin_size_sd = [
        length(collect(skipmissing([row[col] for col in basin_cols]))) > 1 ?
        std(skipmissing([row[col] for col in basin_cols])) : 0.0
        for row in eachrow(combined)
    ]
    
    # Calculate mean and SD for times
    time_cols = [col for col in names(combined) if startswith(string(col), "time")]
    
    combined.time_mean = [
        mean(skipmissing([row[col] for col in time_cols]))
        for row in eachrow(combined)
    ]
    
    combined.time_sd = [
        length(collect(skipmissing([row[col] for col in time_cols]))) > 1 ?
        std(skipmissing([row[col] for col in time_cols])) : 0.0
        for row in eachrow(combined)
    ]
    
    # Select final columns
    select!(
        combined,
        :states,
        :flag,
        :basin_size_mean,
        :basin_size_sd,
        :time_mean,
        :time_sd
    )
    
    return combined
end

"""
    add_frustration_scores(
        df::DataFrame,
        N::Int,
        variables::Vector{String},
        constants::Vector{String},
        node_names::Vector{String},
        interaction_indices::Tuple
    ) -> DataFrame

Add frustration scores for each attractor state using the existing frustration function.
"""
function add_frustration_scores(
    df::DataFrame,
    N::Int,
    variables::Vector{<:AbstractString},
    constants::Vector{<:AbstractString},
    node_names::Vector{<:AbstractString},
    interaction_indices::Tuple
)
    node_order = vcat(variables, constants)
    
    # Create mapping from node_names (in topology) to positions in state vector
    node_positions = [findfirst(==(name), node_order) for name in node_names]
    
    frustrations = Float64[]
    for row in eachrow(df)
        # Parse states string (format: "001" or "001 → 010 → 001")
        states = row.states
        
        # Calculate average frustration across all states in attractor
        frust_values = Float64[]
        for state in states
            
            # Reorder state according to node_positions for topology
            reordered_state = state[node_positions]
            # Call your existing frustration function
            # Assuming signature: frustration(state, interaction_indices; negConv=true)
            frust = frustration(reordered_state, interaction_indices; negConv=true)
            push!(frust_values, frust)
        end
        
        avg_frustration = mean(frust_values)
        push!(frustrations, avg_frustration)
    end
    
    df.frustration = frustrations
    
    return df
end
