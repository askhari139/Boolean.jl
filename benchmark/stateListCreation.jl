nInit = 10000
n_nodes = 25
# benchmarking stateList creation
@btime begin
    stateList = getindex.([rand([-1, 1], nInit) for i in 1:n_nodes], (1:nInit)')
    for i in 1:nInit
        state = stateList[:, i]
    end
end

@btime begin
    stateList = [rand([-1, 1], n_nodes) for i in 1:nInit]
    for i in 1:nInit
        state = stateList[i]
    end
end

# second method is faster by ~2x (3.5ms vs 1.7ms). But can't implement it in shubhamBoolean for partial multilevel 
# benchmarking the difference in indexing
stateMat = getindex.([rand([-1, 1], nInit) for i in 1:n_nodes], (1:nInit)')
stateList = [rand([-1, 1], n_nodes) for i in 1:nInit]
@btime begin
    for i in 1:nInit
        state = stateMat[:, i]
    end
end

@btime begin
    for i in 1:nInit
        state = stateList[i]
    end
end
# indexing is ridiculously faster in the list than matrix. 

@benchmark begin
    stateList = getindex.([rand([-1, 1], nInit) for i in 1:n_nodes], (1:nInit)')
    sL = [stateList[:, i] for i in 1:nInit]
    for i in 1:nInit
        state = sL[i]
    end
end
# 2ms
@benchmark begin
    stateList = [rand([-1, 1], n_nodes) for i in 1:nInit]
    for i in 1:nInit
        state = stateList[i]
    end
end


### check for oeNodes and kdNodes
kdNodes = [2,4,9]
oeNodes = [1,3,5]

@btime begin
    stateList = getindex.([rand([-1, 1], nInit) for i in 1:n_nodes], (1:nInit)')
    if (length(kdNodes) != 0)
        stateList[kdNodes, :] .= stateVec[1]
    end
    if (length(oeNodes) != 0)
        stateList[oeNodes, :] .= stateVec[2]
    end
    sL = [stateList[:, i] for i in 1:nInit]
    for i in 1:nInit
        state = sL[i]
    end
end

@btime begin
    stateList = getindex.([rand([-1, 1], nInit) for i in 1:n_nodes], (1:nInit)')
    if (length(kdNodes) != 0)
        stateList[kdNodes, :] .= stateVec[1]
    end
    if (length(oeNodes) != 0)
        stateList[oeNodes, :] .= stateVec[2]
    end
    for i in 1:nInit
        state = stateList[:, i]
    end
end

@btime begin
    stateList = [rand([-1, 1], n_nodes) for i in 1:nInit]
    for state in stateList
        if !isempty(kdNodes)
            state[kdNodes] .= stateVec[1]  # knockdown: set to -1 or 0
        end
        if !isempty(oeNodes)
            state[oeNodes] .= stateVec[2]  # overexpression: set to +1
        end
    end
    for i in 1:nInit
        state = stateList[i]
    end
end