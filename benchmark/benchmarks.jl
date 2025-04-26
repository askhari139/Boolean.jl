### for sum(Jij*Si*Sj) calculations, we can either perform iterations, or do element-wise multiplications.
### The code below is to test the speed of either methods.
topoFile = "EMT_RACIPE.topo"
uMat, nodes = topo2interaction(topoFile)
nInit = 10000
n_nodes = 25
### Method 1: Iterations
function prodTestIter(uMat, nIter)
    n_nodes = size(uMat,1)
    stateVec = rand([-1,1], n_nodes)
    
    return stateVec
end

# benchmarking zeroConv
stateList = getindex.([rand(stateVec, nInit) for i in 1:n_nodes], (1:nInit)')
state = stateList[:, 1]
@btime join(Boolean.zeroConv($state), "_")
@btime begin
    zc = Boolean.zeroConv($state)
    join(zc, "_")
end

@btime begin
    for i in 1:nInit
        state = stateList[:, i]
        join(Boolean.zeroConv($state), "_")
    end
end

@btime begin
    sL = [stateList[:, i] for i in 1:nInit]
    zc = Boolean.zeroConv.(sL)
    join.(zc, "_")
end
