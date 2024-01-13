## Converts topo file to interaction matrix
function topo2interaction(topoFile::String, type::Int=0)
    df = DataFrame(CSV.File(topoFile))
    dropmissing!(df)
    Nodes = sort(unique(vcat(df.Source, df.Target)))
    n_nodes = length(Nodes)
    update_matrix = zeros(Int64, n_nodes, n_nodes)
    for i in 1:size(df, 1)
        if df[i, 3] == 2
            df[i,3] = -1
        end
        j = findfirst(x->x==df[i,1], Nodes)
        k = findfirst(x->x==df[i,2], Nodes)
        update_matrix[j,k] = Int64(df[i,3])
    end
    if type == 1
        replace!(x -> x == 1 ? 100 : x, update_matrix)
    end

    if type == 2
        replace!(x -> x == -1 ? -100 : x, update_matrix)
    end
    return update_matrix,Nodes
end

## create an exhaustive list of binary states of a given length
function listStates(nStates::Int)
    x = collect(Iterators.product([[-1,1] for i = 1:nStates]...))
    y = []
    for i in x
        j = collect(i)
        push!(y,j)
    end
    return y
end

## table function in R
function freqCalc(x::Array{String,1})
    y = unique(x)
    d = Dict([(i,count(k->k==i,x)) for i in y])
    df = DataFrame(d)
end

## converts binary state (-1s) to decimal
function stateToDec(state::Array{Int,1}, binVec::Array{Int,1})
    y = sum([(i != 1 && i != -1) for i in state])
    if y != 0
        print("Invalid state")
        return
    end
    state = replace(x -> x == -1 ? 0 : x, state)
    if size(state) != size(binVec)
        state = state'
    end
    s = sum(state.*binVec)
    return s
end

## calculates frustration of a state
# function frustration(state::Array{Int,1}, 
#     nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}})
#     frustration = 0
#     nEdges = length(nonZeros[1])
#     for (x,y,v) in zip(nonZeros...)
#         s = state[x]*state[y]*v
#         if s<0
#             frustration = frustration + 1
#         end
#     end
#     frustration = frustration/nEdges
#     return frustration
# end
function frustration(state::Union{Array{Int,1}, Array{Float64,1}}, 
    nonZeros::Union{Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}},
    Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        s = state[x]*state[y]*v
        if s<0
            frustration = frustration + abs(s)
        end
    end
    frustration = frustration/nEdges
    return frustration
end

## 
function frustration0(state::Array{Int,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        p = state[x] == state[y]
        if ((p == true) & (v == 1)) | ((p == false) & (v == -1))
            frustration = frustration + 1
        end
    end
    frustration = frustration/nEdges
    return frustration
end

function getFreq(x)
    y = x/sum(x)
    return y
end

function avg(x)
    m = sum(x)/length(x)
    return m
end

## calculate standard deviation
function SD(x)
    m = avg(x)
    x = [i for i in x]
    # print(x)
    v = sum((x.-m).^2)/length(x)
    s = sqrt(v)
    return s
end

## get frequency for select columns in a dataframe
function dfFreq(state_df::DataFrame, cols::Array{Symbol, 1})
    df = @pipe state_df |>
        groupby(_, cols) |>
        combine(_, nrow => :Count) |>
        transform(_, :Count => getFreq => :frequency) |>
        select(_, push!(cols, :frequency)) |>
        rename(_, :fin => :states)
    return df
end

function dfFreqGen(state_df::DataFrame, cols::Array{Symbol, 1})
    df = @pipe state_df |>
        groupby(_, cols) |>
        combine(_, nrow => :Count) |>
        transform(_, :Count => getFreq => :frequency) |>
        select(_, push!(cols, :frequency))
    return df
end

function dfAvgGen(state_df::DataFrame, cols::Array{Symbol, 1}, meanCol::Array{Symbol, 1})
    df = @pipe state_df |>
        groupby(_, cols) |>
        combine(_, meanCol .=> avg, renamecols = false) |>
        # transform(_, :Count => getFreq => :frequency) |>
        select(_, vcat(cols, meanCol))
    return df
end

## calculate average


## apply a function to each row of a dataframe
function rowWise(df::DataFrame, func::Function)
    v = []
    for i in eachrow(df)
        push!(v, func(Array(i)))
    end
    return v
end

function replaceMissing(x)
    x[ismissing.(x)] .= 0.0
    return x
end

## get mean and SD rowwise of columns containing a keyword in their name
function meanSD(df::DataFrame, keyword::String)
    cols = names(df)
    cols = cols[[occursin(keyword, i) for i in cols]]
    df_new = df[:, cols]
    d = @pipe df |>
        transform(_, cols .=> replaceMissing .=> cols) |>
        transform(_, AsTable(cols) => ByRow(avg) => :Avg) |>
        transform(_, AsTable(cols) => ByRow(SD) => :SD) |>
        select(_, Not(cols))
    return d
end

## convert a state from -1 to 0
function zeroConv(state::Array{Int64,1})
    replace(x -> x == -1 ? 0 : x, sign.(state))
end

function zeroConv(state::Array{Float64,1})
    replace(x -> x == -1 ? 0 : x, Int.(sign.(state)))
end

function getNodes(topoFile::String)
    nodesName = replace(topoFile, ".topo" => "_nodes.txt")
    update_matrix,Nodes = topo2interaction(topoFile)
    io = open(nodesName, "w")
    for i in Nodes
        println(io, i)
    end
    close(io);
end
