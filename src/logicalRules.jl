#=
Author : Kishore Hari
function name : getDegree
Description : This function calculates the in-degree of each node in the network
Inputs :
    topoFile : string; Path to the topology file
    ising : boolean; true if the network is an Ising network, false otherwise
        
Active outputs :
    nodes : list; List of nodes in the network
    intMat : matrix; Interaction matrix of the network
    indegree : list; List of in-degrees of each node

=#

function getDegree(topoFile::String, ising::Bool)
    intMat, nodes = topo2interaction(topoFile)
    if ising
        idMat = Matrix(I, size(intMat)).*0.5
        intMat = intMat + idMat
    end
    # for each column in intMat get the index of non-zero elements
    indegree = [findall(!iszero, intMat[:, j]) for j in 1:size(intMat, 2)]
    return nodes, intMat, indegree
end


#=
function nameCan you make a studio ghibli avatar of a photo? : getUpdateIsing
Description : This function calculates the update function for each node in the network
Inputs :
    topoFile : string; Path to the topology file
Output : 
    indegree : list; List of in-degrees of each node
    updF : list; List of update functions for each node
=#

function getUpdateIsing(topoFile::String)
    nodes, intMat, indegree = getDegree(topoFile, true)
    updF = []
    for (i, node) in enumerate(nodes)
        upd = intMat[:, i]
        upd = upd[upd .!= 0]
        upd = reverse(upd)
        n = length(upd)
        # generate all possible states of length n
        states = Iterators.product(fill([-1, 1], n)...)
        states = reshape(collect.(states), :, 1)
        updVal = [sum(s.*upd) for s in states]
        upF = [ifelse(u < 0, 0, 1) for u in updVal]
        push!(updF, upF)
    end
    return indegree, updF
end

#=
function name : getUpdateNising
Description : This function calculates the update function for each node in the nIsing network
Inputs :
    topoFile : string; Path to the topology file
Output : 
    indegree : list; List of in-degrees of each node
    updF : list; List of update functions for each node
=#


function getUpdateNising(topoFile::String)
    nodes, intMat, indegree = getDegree(topoFile, false)
    updF = []
    for (i, node) in enumerate(nodes)
        upd = intMat[:, i]
        upd = upd[upd .!= 0]
        upd = reverse(upd)
        n = length(upd)
        if n == 0
            push!(updF, [0, 1])
            continue
        end
        # generate all possible states of length n
        states = Iterators.product(fill([0, 1], n)...)
        states = reshape(collect.(states), :, 1)
        updVal = [sum(s.*upd) for s in states]
        upF = [ifelse(u <= 0, 0, 1) for u in updVal]
        push!(updF, upF)
    end
    return indegree, updF
end

#=
function name : getUpdateRand
Description : This function calculates the update function for each node in the random network
Inputs :
    topoFile : string; Path to the topology file

Output :
    indegree : list; List of in-degrees of each node
    updF : list; List of update functions for each node
=#

function getUpdateRand(topoFile::String)
    nodes, intMat, indegree = getDegree(topoFile, false)
    updF = []
    for (i, node) in enumerate(nodes)
        n = length(indegree[i])
        if n < 2
            push!(updF, [0, 1])
            continue
        end
        updVal = rand(0:1, 2^n)
        push!(updF, updVal)
    end
    return indegree, updF
end

#=
function name : getUpdateNCF
Description : This function calculates the update function for each node in the network with NCF
Inputs :
    topoFile : string; Path to the topology file
    ncfOrder : string; Order of NCF
    ord : list; Order of nodes
    
Output :
    indegree : list; List of in-degrees of each node
    updF : list; List of update functions for each node

Notes :
    There are two sources of randomness in the update function:
    1. The order of the nodes
    2. The canalizing variables
=#

function getUpdateNCF(topoFile::String; 
        ncfOrder::String = "random",
        ord::Array{String,1} = String[])
    nodes, intMat, indegree = getDegree(topoFile, false)
    if ncfOrder == "fixed"
        if isempty(ord)
            ord = shuffle(nodes)
        end
        ordNum = [findfirst(x->x==i, nodes) for i in ord]
        ordList = [filter(j -> j in i, ordNum) for i in indegree]
    else
        ordList = [shuffle(i) for i in indegree]
    end
    canVal = [rand(0:1, length(i)) for i in ordList]
    # if the length of canval element is greater than 2, set the last value of the element equal to the second last value
    # canVal = [length(i) > 1 ? vcat(i[1:end-1], i[end-1]) : i for i in canVal]
    updF = []
    for (i, node) in enumerate(nodes)
        print(i)
        upd = intMat[:, i]
        upd1 = upd[upd .!= 0]
        n = length(upd1)
        if n < 2
            push!(updF, [0, 1])
            continue
        end
        states = collect.(Iterators.product(fill([0, 1], n)...))
        states = reshape(states, :, 1)
        states = [reverse(s) for s in states]
        updVal = fill(2, size(states, 1))
        ncfOrd = ordList[i]
        canValNode = canVal[i]
        updat = []
        links = []
        for (j, ncf) in enumerate(ncfOrd)
            cn = canValNode[j]
            link = upd[ncf]
            push!(links, link)
            push!(updat, link < 0 ? 1-cn : cn)
        end
        l = length(updat)
        if l > 1
            if updat[l] != updat[l-1]
                if links[l]*links[l-1] == 1
                    canValNode[l] = canValNode[l-1]
                else
                    canValNode[l] = 1-canValNode[l-1]
                end
                updat[l] = updat[l-1]
            end
        end
        absent = 1-updat[l]
        ncfOrdKey = sortperm(sortperm(ncfOrd))
        # sortperm gives the order of indices that will arrange ncfOrd in a descending order. 
        # Applying sortperm on the result will tell us what is the position of each element of ncfOrd in the states df.
        for (j, ncf) in enumerate(ncfOrdKey)
            id = findall([(updVal[x] == 2) && states[x][ncf] == canValNode[j] for x in 1:size(states, 1)])
            updVal[id] .= updat[j]
        end
        if any(updVal.==2)
            id = findall(updVal .== 2)
            updVal[id] .= absent
        end
        push!(updF, updVal)
    end
    return indegree, updF
end

