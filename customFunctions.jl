### Edge weights

function edgeWeightPert(topoFile::String; nPerts::Int=10000, nInit::Int64=10000, 
    nIter::Int64=1000,
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
        bmodel_reps(topoFile; nInit = nInit, nIter = nIter, mode = mode, stateRep = stateRep, 
        randSim=true, root = string(i), 
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
    root::String="", nLevels = 2, progressMeter::Bool = false, 
    getData::Bool = false)
    update_matrix,Nodes = topo2interaction(topoFile)
    n_nodes = size(update_matrix,1)

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
    Threads.@threads for i in 1:(n_nodes+1)
        if (i == n_nodes + 1)
            turnOffNodes = Int[]
        else
            turnOffNodes = [i]
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
            CSV.write(replace(topoFile, ".topo" => "_singleNodeTurnOff_initFinFlagFreq.csv"), initDf)
        end
        CSV.write(replace(topoFile, ".topo" => "_singleNodeTurnOff_finFlagFreq.csv"), finFlagDf)
    if getData
        if init
            return initDf, finFlagDf
        else
            return finFlagDf
        end
    end
end