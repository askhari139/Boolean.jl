function build_root_name(topoFile; root="", shubham=false, nLevels=2, vaibhav=false,
    turnOffNodes=Int[], turnOffKey="", nNodes=0, oddLevel=false, negativeOdd=false,
    csb=false, timeStep=0.1, stateRep=-1, nonMatrix=false, deleteNodes=Int[],
    kdNodes=Int[], oeNodes=Int[]
)
    rootName = replace(topoFile, ".topo" => "")
    rootName *= (root != "") ? "_$root" : ""
    if shubham
        rootName *= "_shubham_$nLevels"
    end
    if vaibhav
        rootName *= "_turnOff" * (isempty(turnOffNodes) || length(turnOffNodes) == nNodes ? "_All" : (turnOffKey != "" ? "_$turnOffKey" : "_" * join(turnOffNodes, "_")))
    end
    if oddLevel
        rootName *= negativeOdd ? "_oddLevel_negative" : "_oddLevel_positive"
    end
    if csb
        rootName *= "_csb_$timeStep"
    end
    if stateRep == 0
        rootName *= "_nIsing"
    end
    if nonMatrix
        rootName *= "_nonMatrix"
    end
    if !isempty(deleteNodes)
        rootName *= "_del_" * join(deleteNodes, "_")
    end
    if !isempty(kdNodes)
        rootName *= "_kd_" * join(kdNodes, "_")
    end
    if !isempty(oeNodes)
        rootName *= "_oe_" * join(oeNodes, "_")
    end
    return rootName
end





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
function bmodel(topoFile::String; nInit::Int=10_000, nIter::Int=1_000,
    mode::String="Async", stateRep::Int=-1, type::Int=0, randSim::Bool=false,
    randVec::Vector{Float64}=[0.0], shubham::Bool=false, 
    nLevels::Union{Int, Vector{Int}, String}=2,
    vaibhav::Bool=false, csb::Bool=false, timeStep::Float64=0.1,
    discreteState::Bool=true, nonMatrix::Bool=true,
    turnOffNodes::Vector{Int}=Int[], oddLevel::Bool=false, negativeOdd::Bool=false,
    maxNodes::Int=100, kdNodes::Vector{Int}=Int[], oeNodes::Vector{Int}=Int[],
    deleteNodes::Vector{Int}=Int[]
)
    update_matrix, Nodes = topo2interaction(topoFile, type)
    
    if length(Nodes) > maxNodes
        @warn "Too many nodes ($length(Nodes)) > maxNodes ($maxNodes). Exiting."
        return 0, 0, 0
    end

    # Handle deletion of nodes if needed
    if !isempty(deleteNodes)
        keepNodes = setdiff(1:size(update_matrix, 1), deleteNodes)
        update_matrix = update_matrix[keepNodes, keepNodes]
        Nodes = Nodes[keepNodes]
        kdNodes = adjust_indices(setdiff(kdNodes, deleteNodes), deleteNodes)
        oeNodes = adjust_indices(setdiff(oeNodes, deleteNodes), deleteNodes)
    end

    if shubham
        state_df, frust_df = (nLevels == 0) ?
            asyncUpdate(update_matrix, nInit, nIter, stateRep, vaibhav, turnOffNodes, kdNodes, oeNodes) :
            shubhamBoolean(update_matrix, nInit, nIter, nLevels, vaibhav, turnOffNodes, kdNodes, oeNodes)

    elseif csb
        state_df, frust_df = csbUpdate(update_matrix, nInit, nIter; timeStep=timeStep, discreteState=discreteState)

    elseif oddLevel
        state_df, frust_df = oddLevels(update_matrix, nInit, nIter, negativeOdd)

    elseif mode == "Async"
        if nonMatrix
            state_df, frust_df = (stateRep == -1) ? asyncIsingNoFunc(update_matrix, nInit, nIter) :
                                                   asyncNIsingNoFunc(update_matrix, nInit, nIter)
        elseif randSim
            state_df, frust_df = asyncRandUpdate(update_matrix, nInit, nIter, randVec, stateRep)
        else
            state_df, frust_df = asyncUpdate(update_matrix, nInit, nIter, stateRep, vaibhav, turnOffNodes, kdNodes, oeNodes)
        end
    else
        error("Unsupported mode: $mode")
    end

    return state_df, Nodes, frust_df
end

"""
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
"""

function bmodel_reps(topoFile::String; nInit::Int=10_000, nIter::Int=1_000,
    mode::String="Async", stateRep::Int=-1, reps::Int=3, types::Vector{Int}=[0],
    init::Bool=false, randSim::Bool=false, root::String="", randVec::Vector{Float64}=[0.0],
    shubham::Bool=false, nLevels::Union{Int, Vector{Int}, String}=2, vaibhav::Bool=false,
    csb::Bool=false, timeStep::Float64=0.1, discreteState::Bool=false, nonMatrix::Bool=false,
    turnOffNodes::Union{Int, Vector{Int}}=Int[], turnOffKey::String="",
    oddLevel::Bool=false, negativeOdd::Bool=false, write::Bool=true,
    getData::Bool=false, maxNodes::Int=100, kdNodes::Vector{Int}=Int[], 
    oeNodes::Vector{Int}=Int[], deleteNodes::Vector{Int}=Int[]
)
    if isa(turnOffNodes, Int)
        turnOffNodes = [turnOffNodes]
    end

    update_matrix, Nodes = topo2interaction(topoFile)
    nNodes = length(Nodes)

    finFlagFreq_list = []
    initFinFlagFreq_list = []

    for type in types
        single_type_finFlag_list = []
        single_type_initFinFlag_list = []
        frust_list = []

        for rep in 1:reps
            state_df, _, frust_df = bmodel(topoFile, nInit=nInit, nIter=nIter, 
                mode=mode, stateRep=stateRep, type=type, randSim=randSim, randVec=randVec,
                shubham=shubham, nLevels=nLevels, vaibhav=vaibhav, csb=csb,
                timeStep=timeStep, discreteState=discreteState, nonMatrix=nonMatrix,
                turnOffNodes=turnOffNodes, oddLevel=oddLevel, negativeOdd=negativeOdd,
                maxNodes=maxNodes, kdNodes=kdNodes, oeNodes=oeNodes, deleteNodes=deleteNodes)

            if state_df == 0
                error("Too many nodes. Exiting.")
            end

            push!(frust_list, frust_df)
            push!(single_type_finFlag_list, dfFreq(state_df, [:fin, :flag]))
            if init
                push!(single_type_initFinFlag_list, dfFreq(state_df, [:fin, :flag, :init]))
            end
        end

        # Aggregate across reps
        finFlag_df = reduce((x, y) -> outerjoin(x, y, on=[:states, :flag], makeunique=true), single_type_finFlag_list)
        finFlag_df = meanSD(finFlag_df, "frequency")
        frust_df = dfAvgGen(unique(reduce(vcat, frust_list), [:fin, :time]), [:fin, :frust], [:time])
        finFlag_df = outerjoin(finFlag_df, frust_df, on=:states => :fin, makeunique=true)
        rename!(finFlag_df, Dict(:Avg => Symbol("Avg$type"), :SD => Symbol("SD$type"), :frust => Symbol("frust$type")))
        push!(finFlagFreq_list, finFlag_df)

        if init
            initFin_df = reduce((x, y) -> outerjoin(x, y, on=[:init, :states, :flag], makeunique=true), single_type_initFinFlag_list)
            initFin_df = meanSD(initFin_df, "frequency")
            initFin_df = outerjoin(initFin_df, frust_df, on=:states => :fin, makeunique=true)
            rename!(initFin_df, Dict(:Avg => Symbol("Avg$type"), :SD => Symbol("SD$type"), :frust => Symbol("frust$type")))
            push!(initFinFlagFreq_list, initFin_df)
        end
    end

    finFlagFreqFinal_df = reduce((x, y) -> outerjoin(x, y, on=[:states, :flag], makeunique=true), finFlagFreq_list)
    if init
        initFinFlagFreqFinal_df = reduce((x, y) -> outerjoin(x, y, on=[:init, :states, :flag], makeunique=true), initFinFlagFreq_list)
    end

    if vaibhav
        turnOffLabel = isempty(turnOffNodes) ? "All" : join(Nodes[turnOffNodes], "_")
        finFlagFreqFinal_df[!, "turnOffNode"] .= turnOffLabel
        if init
            initFinFlagFreqFinal_df[!, "turnOffNode"] .= turnOffLabel
        end
    end

    if write
        rootName = build_root_name(topoFile; root=root, shubham=shubham, nLevels=nLevels,
            vaibhav=vaibhav, turnOffNodes=turnOffNodes, turnOffKey=turnOffKey, nNodes=nNodes,
            oddLevel=oddLevel, negativeOdd=negativeOdd, csb=csb, timeStep=timeStep,
            stateRep=stateRep, nonMatrix=nonMatrix, deleteNodes=deleteNodes,
            kdNodes=kdNodes, oeNodes=oeNodes)

        finFlagFreqFinal_df = sort(finFlagFreqFinal_df, order(:Avg0, rev=true))
        CSV.write(join([rootName, "_finFlagFreq.csv"]), finFlagFreqFinal_df)

        if init
            initFinFlagFreqFinal_df = sort(initFinFlagFreqFinal_df, order(:Avg0, rev=true))
            CSV.write(join([rootName, "_initFinFlagFreq.csv"]), initFinFlagFreqFinal_df)
        end

        if !randSim
            getNodes(topoFile)
        end
    end

    return getData ? (init ? (finFlagFreqFinal_df, initFinFlagFreqFinal_df) : finFlagFreqFinal_df) : nothing
end


