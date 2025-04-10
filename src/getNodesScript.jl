include(ENV["BMODEL"]*"/bmodel.jl")
using Base.Threads

fileList = readdir()
topoFiles = String[]
for i in fileList
	if endswith(i, "topo")
		push!(topoFiles, i)
	end
end

for topoFile in topoFiles
    getNodes(topoFile)
end