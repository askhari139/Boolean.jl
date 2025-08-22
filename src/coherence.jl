# include("dependencies.jl")
# include("utils.jl")
# include("bmodel.jl")

function charToState(state, nLevels)
    if (nLevels == 0)
        stVec = [-1, 1]
    else
        stVec = getStateVec(nLevels)
    end
    state = parse.(Int, split(state, "_"))
    return [stVec[i + 1] for i in state]
end

function asyncInit(update_matrix2,
    nIter::Int, state::Array{Int,1})
    # state = rand(stateVec, n_nodes) #pick random state
    init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
    flag = 0
    n_nodes = length(state)
    for j in 1:nIter
        s1 = sign.(update_matrix2*state)
        u = rand(1:n_nodes, 1)
        if iszero(j%2) # check after every two steps,hopefully reduce the time
            if s1 == state
                flag = 1
                break
            end
        end
        state[u] = s1[u]
    end
    return state, flag
end

function coherenceIter(update_matrix,
    nIter::Int, state::Union{Vector{Int}, Vector{Float64}},
    pertNodes::Union{Int, Array{Int,1}}, 
    nSim::Int, nLevels::Int)
    init = copy(state)
    state[pertNodes] = -1*sign.(state[pertNodes])
    stateList = [state for i in 1:nSim]
    stateList = vcat([init], stateList)
    if (nLevels == 0)
        state_df, frust_df = asyncUpdate(update_matrix, 0, nIter, -1, false, Int[], Int[], Int[]; stateList = stateList)
    else 
        state_df, frust_df = shubhamBoolean(update_matrix, 0, nIter, nLevels, false, Int[], Int[], Int[]; stateList = stateList)
    end
    stInit = state_df[1,1]
    state_df = state_df[2:end,:]
    # frust_df = frust_df[2:end, :]
    state_df = state_df |>
        x -> rename!(x, :init => :pert) |>
        x -> insertcols!(x, :init => [stInit for i in 1:size(state_df, 1)]) |>
        x -> select!(x, [:init, :pert, :fin, :flag]) |>
        x -> leftjoin(x, frust_df, on = :fin)
    return state_df
end

hamming_str(s1, s2) = sum(parse.(Int, split(s1, "_")) .!= parse.(Int, split(s2, "_")))/length(split(s1, "_"))

function coherence(topoFile::String; nIter::Int=1000, 
    nPert::Int=1, nInit::Int=100, nSim::Int=10, nLevels::Int=1, 
    rowz::Vector{Int} = Vector{Int}(), returnDf = false)
    update_matrix,Nodes = topo2interaction(topoFile)
    n_nodes = length(Nodes)
    if (nLevels == 0)
        key = "_finFlagFreq.csv"
    else
        key = "_shubham_"*string(nLevels)*"_finFlagFreq.csv"
    end
    fileName = replace(topoFile, ".topo" => key)
    states = CSV.read(fileName, DataFrame) 
    states = states[:, :states]
    rowz = filter(x -> x < length(states), rowz)
    states = states[rowz]
    states = [charToState(st, nLevels) for st in states]
    b = binomial(n_nodes, nPert)
    if b < 2*nInit
        nInit = min(nInit, b)
        choice = collect(combinations(1:n_nodes, nPert))[rand(1:b, nInit)]
        choice = vcat((x' for x in choice)...)
    else
        choice = reshape(rand(1:n_nodes, nPert*nInit), nInit, nPert)
    end
    collectionLarge = DataFrame(init=String[], pert = String[], fin=String[], flag=Int[], frust=Float64[], time = Float64[])
    for st in states
        for i in 1:nInit
            pertNodes = choice[i,:]
            collectionLarge = vcat(collectionLarge, 
                coherenceIter(update_matrix, nIter, st, pertNodes, nSim, nLevels))
        end
    end

    collectionLarge.nPert = fill(nPert, size(collectionLarge, 1))
    collectionLarge.returned = Int.(collectionLarge.init .== collectionLarge.fin)
    collectionLarge.hamming = hamming_str.(collectionLarge.init, collectionLarge.fin)
    nm = join(["_nPert_", string(nPert), 
        "_nLevs_", string(nLevels) , "_coherence.csv"])
    rootName = replace(topoFile, ".topo" => nm)
    CSV.write(rootName, collectionLarge)
    if returnDf
        return collectionLarge
    end
end


function coherenceAllNode(topoFile::String; nIter::Int=1000, 
    nInit::Int=100, nSim::Int=10, nLevels::Int=1, 
    rowz::Vector{Int} = Vector{Int}())
    update_matrix,Nodes = topo2interaction(topoFile)
    n_nodes = length(Nodes)
    dfList = [coherence(topoFile; nIter = nIter, nPert = i, 
        nInit = nInit, nSim = nSim,
            nLevels = nLevels, rowz = rowz, returnDf = true) for i in 1:n_nodes]
    df = vcat(dfList...)
    df1 = df |> 
        x -> begin
        # Convert to Int/float and replace missing with 0
            x.returned = Int.(coalesce.(x.returned, 0))
            x.hamming = float.(coalesce.(x.hamming, 0))
            x
        end |>
        x -> groupby(x, [:init, :nPert]) |>
        g -> combine(g,
                    :hamming => avg => :hamming_mean,
                    :hamming => SD  => :hamming_sd,
                    :returned => avg => :returned_mean,
                    :returned => SD  => :returned_sd)
    rootName = replace(topoFile, ".topo" => "_nLevels_" * string(nLevels) * "_coherence.csv")
    CSV.write(rootName, df1)
end

# function coherence(topoFile::String; nIter::Int=1000, simulate::Bool=false,
#     bmodelArgs::Dict=Dict(), nPert::Int=1, randPert::Bool=false, 
#     nInit::Int=100, nSim::Int=10, nLevels::Int=1, rowz::Vector{Int} = Vector{Int}())
#     if simulate
#         x = bmodel_reps(topoFile, bmodelArgs...)
#     end
#     update_matrix,Nodes = topo2interaction(topoFile)
#     n_nodes = length(Nodes)
#     update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
#     update_matrix2 = sparse(update_matrix2')
#     fileName = replace(topoFile, ".topo" => "_finFlagFreq.csv")
#     states = CSV.read(fileName, DataFrame) 
#     states = filter(:flag => x -> x== 1, states)
#     if typeof(states.Avg0) == Vector{String}
#         states = filter(:Avg0 => x -> x!="NA", states)
#     else
#         states = filter(:Avg0 => x -> !ismissing(x), states)
#     end
#     states = states[:, :states]
#     nNodes = length(states[1]) - 2
#     if nPert > 1
#         randPert = true
#     end
#     choice = reshape(rand(1:nNodes, nPert*nInit), nInit, nPert)
#     if !randPert
#         choice = collect(1:nNodes)
#         nInit = nNodes
#     end
#     collectionLarge = DataFrame(init=[], fin=[], flag=[], frequency=Float64[])
#     for st in states
#         st = st[2:(nNodes+1)]
#         st = [x=='0' ? -1 : 1 for x in st]
#         for i in 1:nInit
#             pertNodes = choice[i,:]
#             collectionLarge = vcat(collectionLarge, coherenceIter(update_matrix2, nIter, st, pertNodes, nSim))
#         end
#     end
#     nm = "singleNode"
#     if !randPert
#         nm = "allNode"
#     end
#     if nPert > 1
#         nm = "multiNode"
#     end
#     nm = join([nm, "_coherence.csv"])
#     rootName = replace(topoFile, ".topo" => nm)
#     CSV.write(rootName, collectionLarge)
# end

