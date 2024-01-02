### Edge weights

function edgeWeightPert(topoFile::String; nPerts::Int=10000, nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, reps::Int = 3, csv::Bool=false, 
    types::Array{Int, 1} = [0,1,2],init::Bool=false, randSim::Bool=true)
    updMat, nodes = topo2interaction(topoFile)
    nZ = length(findall(updMat.!=0))
    nRand = nZ*nPerts
    rands = reshape(rand(nRand), nPerts, nZ)
    randFold = replace(topoFile, ".topo" => "_rand")
    d1 = pwd()
    mkpath(randFold)
    cpPath = join([randFold, "/", topoFile])
    cp(topoFile, cpPath, force = true)
    cd(randFold)
    Threads.@threads for i in 1:nPerts
        # println(string(i))
        bmodel_reps(topoFile; nInit = nInit, nIter = nIter, mode = mode, stateRep = stateRep, randSim=true, root = string(i), 
        randVec = rands[i,:], types = types)
    end
    nodesName = join([replace(topoFile, ".topo" => ""), "_nodes.txt"])
        update_matrix,Nodes = topo2interaction(topoFile)
        io = open(nodesName, "w")
        for i in Nodes
            println(io, i)
        end
    close(io);
    cd(d1)

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
function singleNodeTurnOff(topoFile::String; nInit::Int64=10000, nIter::Int64=1000,
    mode::String="Async", stateRep::Int64=-1, reps::Int = 3, init::Bool=false, 
    root::String="", nLevels = 2)
    update_matrix,Nodes = topo2interaction(topoFile)
    n_nodes = size(update_matrix,1)

    ### Normal simulation
    dfs = bmodel_reps(topoFile; nInit = nInit, nIter = nIter, 
            mode = mode, stateRep = stateRep, reps = reps, init = init, nLevels = nLevels,
            root = root, shubham = true, vaibhav = false, write = false)
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
    for i in 1:(n_nodes+1)
        if (i == n_nodes + 1)
            turnOffNodes = Int[]
        else
            turnOffNodes = [i]
        end
        dfs = bmodel_reps(topoFile; nInit = nInit, nIter = nIter, 
            mode = mode, stateRep = stateRep, reps = reps, init = init, nLevels = nLevels,
            root = root, shubham = true, vaibhav = false, turnOffNodes = turnOffNodes, write = false)
        if init
            finFlag, init = dfs
            push!(initList, init)
        else
            finFlag = dfs
        end
        push!(finFlagList, finFlag)
        if init
            initDf = vcat(initList...)
        end
        finFlagDf = vcat(finFlagList...)
        if init
            CSV.write(replace(topoFile, ".topo" => "_singleNodeTurnOff_init.csv"), initDf)
        end
        CSV.write(replace(topoFile, ".topo" => "_singleNodeTurnOff_finFlag.csv"), finFlagDf)
    end
end