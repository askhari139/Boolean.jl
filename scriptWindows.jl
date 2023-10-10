include("F:/Github/Projects/Ongoing/Boolean_project/shubhamBoolean/Bmodel/bmodel.jl")
using Base.Threads
topoFiles = map(x->string(x), ARGS)
println(Threads.nthreads())

Threads.@threads for topoFile in topoFiles
 	y1 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = -1, shubham=true, nLevels=2)
	y2 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = -1)
	y3 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = 0)
 	y1 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = -1, shubham=true, nLevels=1, vaibhav=true)
 	println(topoFile, " - ", y1, " seconds, ", y2, " seconds.")
end
