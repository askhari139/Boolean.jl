include("/mnt/d/Github/Projects/Ongoing/Relative_boolean/bm_julia_package/bmodel.jl")
using Base.Threads
topoFiles = map(x->string(x), ARGS)
println(Threads.nthreads())

Threads.@threads for topoFile in topoFiles
 	y1 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = -1)
 	#y2 = @elapsed x = bmodel_reps(topoFile; nInit = 100000, nIter = 1000, mode = "Async", stateRep = 0)
 	println(topoFile, " - ", y1, " seconds.")
end
