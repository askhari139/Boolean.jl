# include("/path/to/Bmodel/bmodel.jl")
using Base.Threads
include("bmodel.jl")

minWt = 0.0
maxWt = 0.01
nPert = 10 # Number of samples of edge weights

fileList = readdir()
topoFiles = String[]
for i in fileList
	if endswith(i, "topo")
		push!(topoFiles, i)
	end
end

println(Threads.nthreads())

for topoFile in [topoFiles[3]]
	y2 = @elapsed xVe = edgeWeightPert(topoFile; nPerts = nPert,
			minWeight = minWt, maxWeight = maxWt)
end
