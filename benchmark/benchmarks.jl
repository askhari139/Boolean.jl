### for sum(Jij*Si*Sj) calculations, we can either perform iterations, or do element-wise multiplications.
### The code below is to test the speed of either methods.
topoFile = "EMT_RACIPE.topo"
uMat, nodes = topo2interaction(topoFile)
### Method 1: Iterations
function prodTestIter(uMat, nIter)
    n_nodes = size(uMat,1)
    stateVec = rand([-1,1], n_nodes)
    
    return stateVec
end