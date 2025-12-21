# include("/path/to/Bmodel/bmodel.jl")
using Boolean

weightList = [(0,0.25), (0, 0.5), (0, 0.75), (0, 1), (0.25, 1), (0.5, 1), (0.75, 1)]
nPert = 1000 # Number of samples of edge weights

noiseVals = [0.001, 0.005, 0.01, 0.02, 0.05]


fileList = readdir()
topoFiles = String[]
for i in fileList
	if endswith(i, "topo")
		push!(topoFiles, i)
	end
end

println(Threads.nthreads())
for noise in noiseVals
	for topoFile in topoFiles
		y2 = @elapsed xVe = contWeightPert(topoFile; nInit=1000,
    nIter=100000, mode="Async", stateRep=-1, 
    noise=noise, steadyStates=true, topN=10)
		println(y2)
	end
end
