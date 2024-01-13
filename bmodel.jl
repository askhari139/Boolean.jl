include("dependencies.jl")
include("utils.jl")
include("async_update.jl")
# include("shubhamBoolean.jl")
include("multiLevel_shubham.jl")
include("CSB.jl")
include("async_non_matrix.jl")
include("customFunctions.jl")
include("oddLevel.jl")

#=
Author : Kishore Hari
function name : bmodel
Inputs :
    topoFile : string; path to the topo file for the network
    nInit : integer; number of initial conditions to use for simulation
    nIter : integer; maximum number of time steps for simulation
    mode : string; ['Async', 'Sync']; mode of update
    stateRep : integer; [-1, 0]; whether to use -1 or 0 to represent low-expression
    rep : integer; replicate number. 
    csv : boolean; whether to write the output table into csv
    type : integer; [0, 1, 2]; the form of update rule to be used. 0 is equal weightage, 1 is activation dominant and 2 is inhibition dominant
Active outputs :
    state_df : DataFrame; 3 columns - [init, fin, flag]:
        init : string; Initial state
        fin : string; Final state
        flag : integer; 1 - fin is a steady state, 0 - not
    Nodes : string array; List of nodes with the order in which state is listed
Passive outputs :
    The function writes states_df to a csv file if csv = true in input
=#
function bmodel(topoFile::String; nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, type::Int=0, randSim::Bool = false,
    randVec::Array{Float64, 1}=[0.0], shubham = false, discrete = true, nLevels = 2,
    vaibhav::Bool = false, csb::Bool = false, timeStep::Float64 = 0.1,
    discreteState::Bool = true, nonMatrix::Bool = true,
    turnOffNodes::Array{Int,1} = Int[], oddLevel::Bool = false,
    negativeOdd::Bool = false)
    update_matrix,Nodes = topo2interaction(topoFile, type)
    if shubham == true
        state_df, frust_df = shubhamBoolean(update_matrix, nInit, nIter, discrete; nLevels = nLevels, vaibhav = vaibhav, turnOffNodes = turnOffNodes)
    elseif csb == true
        state_df, frust_df = csbUpdate(update_matrix, nInit, nIter; timeStep = timeStep, discreteState = discreteState)
    elseif oddLevel == true
        state_df, frust_df = oddLevels(update_matrix, nInit, nIter, negativeOdd)
    elseif mode == "Async"
        if nonMatrix
            if stateRep == -1
                state_df, frust_df = asyncIsingNoFunc(update_matrix, nInit, nIter)
            else
                state_df, frust_df = asyncNIsingNoFunc(update_matrix, nInit, nIter)
            end
        else
            if stateRep == -1
                if randSim
                    state_df, frust_df = asyncRandUpdate(update_matrix, nInit, nIter, randVec)
                else 
                    state_df, frust_df = asyncUpdate(update_matrix, nInit, nIter)
                end
            else
                state_df, frust_df = asyncUpdate2(update_matrix, nInit, nIter)
            end
        end
    else
        print("Method under construction.")
    end
    # file_name = join([replace(topoFile, ".topo" => "_"), repl])
    # CSV.write(join(name,"_bmRes.csv"]), state_df)
    return state_df,Nodes, frust_df
end

#=
Author : Kishore Hari
function name : bmodel_reps
Inputs :
    topoFile : string; path to the topo file for the network
    nInit : integer; number of initial conditions to use for simulation
    nIter : integer; maximum number of time steps for simulation
    mode : string; ['Async', 'Sync']; mode of update
    stateRep : integer; whether to use -1 or 0 to represent low-expression
    reps : integer; number of replicates 
    csv : boolean; whether to write the output table of the function bmodel into csv
    types : integer array; subset of [0, 1, 2]; the forms of update rule to be used. 0 is equal weightage, 1 is activation dominant and 2 is inhibition dominant
    init : bool; Whether or not to include the initial conditions in the output. 
        If checked true, output will contain the frequency of all unique initial condition-steady state pair. For asynchronous boolean, this increases the table size by a lot
Active outputs :
    finFreqFinal_df : DataFrame; [states, Avg0, SD0, frust0,...]:
        states : string; Steady states
        Avg0 : float; Mean frequency (from 'reps' replicates) of the steady states obtained using rule 0
        SD0 : float; Standard deviation of the frequency
        frust0 : float; Frustration of the steady state
        Similarly we have Avg1, SD1, frust1, Avg2, SD2 and frust 2 depending upon the types argument for the function
    finFlagFreqFinal_df : DataFrame; [states, flag, Avg0,...]
        flag : integer; 1 - states is a steady state, 0 - not
    initFinFlagFreqFinal_df : DataFrame; [init, states, flag, Avg0,...]
        init : string; initial state
Passive outputs :
    The function writes finFreqFinal_df, finFlagFreqFinal_df, initFinFlagFreqFinal_df to files.
=#

function bmodel_reps(topoFile::String; nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, reps::Int = 3,
    types::Array{Int, 1} = [0],init::Bool=false, randSim::Bool=false, root::String="", 
    randVec::Array{Float64,1}=[0.0], shubham = false, discrete = false, nLevels = 2,
    vaibhav::Bool = false, csb::Bool = false, timeStep::Float64 = 0.1,
    discreteState::Bool = false, nonMatrix::Bool = false,
    turnOffNodes::Union{Int64, Array{Int,1}} = Int[], 
    oddLevel::Bool = false, negativeOdd::Bool = false,
    write::Bool = true, getData::Bool = false)
    update_matrix,Nodes = topo2interaction(topoFile)
    nNodes = length(Nodes)
    finFlagFreqFinal_df_list_list = []
    initFinFlagFreqFinal_df_list_list = []
    frust_df_list = []

    if (typeof(turnOffNodes) == Int64)
        turnOffNodes = [turnOffNodes]
    end

    for type in types

        finFlagFreqFinal_df_list = []
        initFinFlagFreqFinal_df_list = []
        frust_df_list = []

        for rep in 1:reps
            states_df, Nodes, frust_df = bmodel(topoFile, nInit = nInit, 
                nIter = nIter, mode = mode, stateRep = stateRep, type = type, 
                randSim = randSim, randVec = randVec, shubham = shubham, 
                discrete = discrete, nLevels = nLevels, vaibhav = vaibhav,
                csb = csb, timeStep = timeStep, discreteState = discreteState,
                nonMatrix = nonMatrix, turnOffNodes = turnOffNodes,
                oddLevel = oddLevel, negativeOdd = negativeOdd)
            # state_df = dropmissing(state_df, disallowmissing = true)
            push!(frust_df_list, frust_df)
            # Frequnecy table 
            #final states with flag
            finFlagFreq_df = dfFreq(states_df, [:fin, :flag])

            # all counts
            if init
                initFinFlagFreq_df = dfFreq(states_df, [:fin, :flag, :init])
                push!(initFinFlagFreqFinal_df_list, initFinFlagFreq_df)
            end
            push!(finFlagFreqFinal_df_list, finFlagFreq_df)
        end

        # println(typeof(finFreqFinal_df))
        finFlagFreqFinal_df = finFlagFreqFinal_df_list[1]
        if init
            initFinFlagFreqFinal_df = initFinFlagFreqFinal_df_list[1]
        end
        for i in 2:reps
            finFlagFreqFinal_df = outerjoin(finFlagFreqFinal_df, 
                finFlagFreqFinal_df_list[i], 
                on = [:states, :flag], makeunique=true)
            if init
                initFinFlagFreqFinal_df = outerjoin(initFinFlagFreqFinal_df, 
                    initFinFlagFreqFinal_df_list[i],
                    on = [:init, :states, :flag], makeunique = true)
            end
        end

        frust_df = reduce(vcat, frust_df_list)
        # for i in frust_df_list
        #     frust_df = vcat(frust_df, i)
        # end
        frust_df = unique(frust_df, [:fin, :time])
        frust_df = dfAvgGen(frust_df, [:fin, :frust], [:time])

        finFlagFreqFinal_df = meanSD(finFlagFreqFinal_df, "frequency")
        finFlagFreqFinal_df = outerjoin(finFlagFreqFinal_df, frust_df, 
            on = :states => :fin, makeunique =true)
        finFlagFreqFinal_df = rename(finFlagFreqFinal_df, 
            Dict(:Avg => Symbol(join(["Avg", type])), 
                :SD => Symbol(join(["SD", type])),
                :frust => Symbol(join(["frust", type]))))
        push!(finFlagFreqFinal_df_list_list, finFlagFreqFinal_df)


        if init
            initFinFlagFreqFinal_df = meanSD(initFinFlagFreqFinal_df, "frequency")
            initFinFlagFreqFinal_df = outerjoin(initFinFlagFreqFinal_df, frust_df, 
                on = :states => :fin, makeunique =true)
            initFinFlagFreqFinal_df = rename(initFinFlagFreqFinal_df, 
                Dict(:Avg => Symbol(join(["Avg", type])), 
                :SD => Symbol(join(["SD", type])),
                :frust => Symbol(join(["frust", type]))))
            push!(initFinFlagFreqFinal_df_list_list, initFinFlagFreqFinal_df)

        end
    end
        # println(typeof(finFreqFinal_df))
    finFlagFreqFinal_df = finFlagFreqFinal_df_list_list[1]
    if init
        initFinFlagFreqFinal_df = initFinFlagFreqFinal_df_list_list[1]
    end
    n = length(types)
    if n > 1
        for i in 2:n
            finFlagFreqFinal_df = outerjoin(finFlagFreqFinal_df, 
                finFlagFreqFinal_df_list_list[i], 
                on = [:states, :flag], makeunique=true)
            if init
                initFinFlagFreqFinal_df = outerjoin(initFinFlagFreqFinal_df, 
                    initFinFlagFreqFinal_df_list_list[i],
                    on = [:init, :states, :flag], makeunique = true)
            end
        end
    end

    if !randSim
        nodesName = replace(topoFile, ".topo" => "_nodes.txt")
        update_matrix,Nodes = topo2interaction(topoFile)
        io = open(nodesName, "w")
        for i in Nodes
            println(io, i)
        end
        close(io);
    end

    if vaibhav
        if (length(turnOffNodes) == 0)
            turnOffLabel = "All"
        else
            turnOffLabel = join(Nodes[turnOffNodes], "_")
        end
        finFlagFreqFinal_df[!, "turnOffNode"] .= turnOffLabel
        if init
            initFinFlagFreqFinal_df[!, "turnOffNode"] .= turnOffLabel
        end
    end

    if write
        rootName = replace(topoFile, ".topo" => "")
        if root !=""
            rootName = join([rootName, "_",root])
        end
        if shubham
            rootName = join([rootName, "_shubham_", nLevels])
            if vaibhav
                rootName = join([rootName, "_turnOff"])
                nTurnOff = length(turnOffNodes)
                if (nTurnOff == nNodes || nTurnOff == 0)
                    rootName = join([rootName, "_All"])
                else
                    tList = join(turnOffNodes, "_")
                    rootName = join([rootName, "_", tList])
                end
            end
        end
        if oddLevel
            rootName = join([rootName, "_oddLevel"])
            if negativeOdd
                rootName = join([rootName, "_negative"])
            else
                rootName = join([rootName, "_positive"])
            end
        end
        if csb
            rootName = join([rootName, "_csb_", timeStep])
        end
        # println(rootName)
        if stateRep == 0
            rootName = join([rootName, "_nIsing"])
        end
        if nonMatrix
            rootName = join([rootName, "_nonMatrix"])
        end
        finFlagFreqName = join([rootName, "_finFlagFreq.csv"])

        CSV.write(finFlagFreqName, finFlagFreqFinal_df)


        if init
            initFinFlagFreqName = join([rootName, "_initFinFlagFreq.csv"])
            CSV.write(initFinFlagFreqName, initFinFlagFreqFinal_df)
        end
    end

    if getData
        if init
            return(finFlagFreqFinal_df, 
                initFinFlagFreqFinal_df)
        else
            return(finFlagFreqFinal_df)
        end
    end
end





