# include("/path/to/Bmodel/bmodel.jl")
using Base.Threads

minWt = 0.0
maxWt = 1.0
nPert = 10 # Number of samples of edge weights

fileList = readdir()
topoFiles = String[]
for i in fileList
	if endswith(i, "topo")
		push!(topoFiles, i)
	end
end

println(Threads.nthreads())

Threads.@threads for topoFile in topoFiles
 	y1 = @elapsed x = edgeWeightPert(topoFile; nPerts = nPert, nInit = 100000,
			minWeight = minWt, maxWeight = maxWt)
end
#for topoFile in topoFiles
#	y3 = @elapsed x = edgeWeightPert(topoFile; nPerts = 100, nInit = 10000, types = [0])
#	println(topoFile, " - ", y3, " seconds.")
#end

