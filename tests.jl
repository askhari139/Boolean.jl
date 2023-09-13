include("bmodel.jl")
topoFile = "EMT_RACIPE.topo"
# y1 = @elapsed x = bmodel_reps(topoFile, nInit = 1000, nIter = 1000, mode = "Async", stateRep = 0, randSim=false, shubham = false, nonMatrix = true)
# y2 = @elapsed x = bmodel_reps(topoFile, nInit = 1000, nIter = 1000, mode = "Async", stateRep = 0, randSim=false, shubham = true, nonMatrix = false, nLevels = 1)
y3 = @elapsed x = bmodel_reps(topoFile, nInit = 1000000, nIter = 1000, mode = "Async", stateRep = 0, randSim=false, shubham = false, nonMatrix = false, init = true)

# println(topoFile, "\nnon matrix: ", y1, " seconds\nshubham: ", y2, " seconds\nIsing: ", y3, " seconds.")
