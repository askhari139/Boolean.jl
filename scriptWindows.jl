include("/mnt/external4tb/Kishore/shubhamBoolean/Bmodel/bmodel.jl")
using Base.Threads
topoFiles = map(x->string(x), ARGS)
println(Threads.nthreads())

Threads.@threads for topoFile in topoFiles
 	y1 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = -1, randSim = false, shubham=true)
 	#y2 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = 0, randSim = false)
 	println(topoFile, " - ", y1, " seconds.")
end
