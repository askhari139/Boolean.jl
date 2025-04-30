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
    p = Progress(nPerts)
    Threads.@threads for i in 1:nPerts
        # println(string(i))
        bmodel_reps(topoFile; nInit = nInit, nIter = nIter, mode = mode,
        stateRep = stateRep, randSim=true, root = string(i), 
        randVec = rands[i,:], types = types, reps = reps)
    end
    finish!(p)
    nodesName = join([replace(topoFile, ".topo" => ""), "_nodes.txt"])
        update_matrix,Nodes = topo2interaction(topoFile)
        io = open(nodesName, "w")
        for i in Nodes
            println(io, i)
        end
    close(io);
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
    noise::Float64=0.01, randVec = [0.0], steadyStates::Bool=true)
    updMat, nodes = topo2interaction(topoFile)
    stateMat, states = asyncRandCont(updMat, nInit, nIter, stateRep; 
        weightFunc = defaultWeightsFunction(noise), randVec = randVec, 
        steadyStates=steadyStates)
    MRT = []
    MRTsd = []
    switchingEvents = []
    for i in eachindex(states)
        m = [sum(stateMat[j, :].== i)/nIter for j in 1:size(stateMat, 1)]
        push!(MRT, avg(m))
        push!(MRTsd, SD(m))
    end
    for i in 1:size(stateMat, 1)
        s = stateMat[i, :]
        counter = 0
        for t in eachindex(s)
            if t == 1
                continue
            end
            if s[t] != s[t-1]
                counter = counter + 1
            end
        end
        push!(switchingEvents, counter)
    end
    # write stateMat to a file with the name topoFile_contWeightPert.dat
    stateMat = DataFrame(stateMat, :auto)
    CSV.write(replace(topoFile, ".topo" => "_contWeightPert.csv"), stateMat)
    # write states to a text file with the name topoFile_contWeightPert_states.txt
    stateDf = DataFrame(states = states, MRT = MRT, MRTsd = MRTsd, 
        ID = 1:length(states))
    CSV.write(replace(topoFile, ".topo" => "_contWeightPertMRT.csv"), stateDf)
    return 
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
function simulate_network_logical(rules::String; 
    folder::String=".",
    reps::Int=3, 
    # text_to_BN args
    max_N::Int=10000,
    # num_of_steady_states_asynchronous args
    nIter::Int=1000, 
    nInit::Int=10000, 
    seed::Union{Int,Nothing}=-1
)
    
    dataFile = replace(rules, ".txt" => ".json")
    if (!isfile(joinpath(folder, dataFile)))
        println("File $dataFile not found in folder $folder")
        l = pyCheck()
        if !l
            println("Python modules not found. Please install them.")
            return
        else
            logicRules.text_to_BN(rules)
            logicRules.boolean_to_topo(rules)
            if (!isfile(joinpath(folder, dataFile)))
                println("Python script failed to create $dataFile.")
                return
            end
        end
    end
    F, I, N, degree, variables, constants = load_bn_from_json(joinpath(folder, dataFile))
    getFrust = true
    topoFile = replace(rules, ".txt" => ".topo")
    if !isfile(joinpath(folder, topoFile))
        println("File $topoFile not found in folder $folder")
        getFrust = false
    end

    if getFrust
        update_matrix, Nodes = topo2interaction(topoFile)
        nzID = findnz(sparse(update_matrix))
    end

    all_runs = Dict{String, Tuple{Int, Float64}}()
    total_points = 0
    dictF = Dict{Tuple{Int, Int}, Int}()
    ssDfList = []
    for _ in 1:reps
        y = @elapsed ssDf, dictF = num_of_steady_states_asynchronous(F, I, N; 
            nsim=nInit, search_depth=nIter, DEBUG = false, SEED=seed, dictF=dictF)
        ssDfList = vcat(ssDfList, ssDf)
    end
    ssDf = reduce((x, y) -> outerjoin(x, y, on=[:steady_state, :flag], makeunique=true), ssDfList)
    states = ssDf[!, :steady_state]
    states = [dec2binvec(i, N) for i in states]
    nodesLog = vcat(variables, constants)
    p = [findfirst(==(x), nodesLog) for x in Nodes]
    if getFrust
        frust = [frustration(state[p], nzID; negConv = true) for state in states]
    else
        frust = [NaN for _ in states]
    end
    states = [join(state[p], "_") for state in states]
    ssDf[!, :states] .= states
    ssDf[!, :frust0] .= frust
    ssDf = meanSD(ssDf, "basin_size"; avgKey=:Avg0, sdKey=:SD0)
    ssDf = meanSD(ssDf, "time"; avgKey=:time, sdKey=:SDTime)
    ssDf = select(ssDf, [:states, :flag, :Avg0, :SD0, :frust0, :time])
    ssDf = sort(ssDf, order(:Avg0, rev=true))
    outFile = replace(rules, ".txt" => "_finFlagFreq.csv")
    CSV.write(joinpath(folder, outFile), ssDf)
end

