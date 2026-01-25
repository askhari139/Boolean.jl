## Converts topo file to interaction matrix
function topo2interaction(topoFile::String, type::Int=0)
    df = DataFrame(CSV.File(topoFile))
    dropmissing!(df)
    Nodes = sort(unique(vcat(df[:,1], df[:,2])))
    n_nodes = length(Nodes)
    types = sort(unique(df[:,3]))
    if (length(types) > 2 && typeof(df[1,3]) == Int) || typeof(df[1,3]) == Float64
        df[:, 3] = float(df[:, 3])
    else
        df[:,3] = replace(x -> x == 2 ? -1 : x, df[:,3])
    end
    update_matrix = zeros(typeof(df[1,3]), n_nodes, n_nodes)
    for row in eachrow(df)
        j = findfirst(x->x==row[1], Nodes)
        k = findfirst(x->x==row[2], Nodes)
        update_matrix[j,k] = row[3]
    end
    if type == 1
        replace!(x -> x >0 ? 100*x : x, update_matrix)
    end

    if type == 2
        replace!(x -> x <0 ? 100*x : x, update_matrix)
    end
    return update_matrix,String.(Nodes)
end

function topo2FIN(topo_file::String)
    up_mat, Nodes = topo2interaction(topo_file)
    I = [sum(x .!=0) for x in eachcol(up_mat)]
    N = length(Nodes)
    return up_mat, I, N
end


## create an exhaustive list of binary states of a given length
function listStates(nStates::Int, stateVec)
    x = collect(Iterators.product([stateVec for i = 1:nStates]...))
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

"""
    bin2dec(bits::Vector{Int}) -> Int

Converts a binary vector (e.g., [1, 0, 1]) to its decimal representation.
"""
function bin2dec(bits::Vector{Int})
    return foldl((acc, b) -> 2 * acc + b, bits)
end

## converts binary state (-1s) to decimal
function stateToDec(state::Array{Int,1}, binVec::Array{Int,1})
    y = sum([(i != 1 && i != -1) for i in state])
    if y != 0
        print("Invalid state")
        return
    end
    state = replace(x -> x == -1 ? 0 : x, state)
    return bin2dec(state)
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
    Tuple{Array{Int64,1},Array{Int64,1},Array{Int64,1}}};negConv = false)
    frustration = 0
    nEdges = length(nonZeros[1])
    if (negConv)
        state = replace(x -> x == 0 ? -1 : x, state)
    end
    for (x,y,v) in zip(nonZeros...)
        s = state[x]*state[y]*v
        if s<0
            frustration = frustration + abs(s)
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
        DataFrames.groupby(_, cols) |>
        combine(_, nrow => :Count) |>
        transform(_, :Count => getFreq => :frequency) |>
        select(_, push!(cols, :frequency)) |>
        rename(_, :fin => :states)
    return df
end

function dfFreqGen(state_df::DataFrame, cols::Array{Symbol, 1})
    df = @pipe state_df |>
        DataFrames.groupby(_, cols) |>
        combine(_, nrow => :Count) |>
        transform(_, :Count => getFreq => :frequency) |>
        select(_, push!(cols, :frequency))
    return df
end

function dfAvgGen(state_df::DataFrame, cols::Array{Symbol, 1}, meanCol::Array{Symbol, 1})
    df = @pipe state_df |>
        DataFrames.groupby(_, cols) |>
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
function meanSD(df::DataFrame, keyword::String; avgKey::Symbol = :Avg, sdKey::Symbol = :SD)
    cols = names(df)
    cols = cols[[occursin(keyword, i) for i in cols]]
    df_new = df[:, cols]
    d = @pipe df |>
        transform(_, cols .=> replaceMissing .=> cols) |>
        transform(_, AsTable(cols) => ByRow(avg) => avgKey) |>
        transform(_, AsTable(cols) => ByRow(SD) => sdKey) |>
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

function signVec(state)
    sign.(state)
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


function getPeripherals(topoFile::String; repeat::Bool=false)
    update_matrix,Nodes = topo2interaction(topoFile)
    # change update_matrix to _positive
    update_matrix = abs.(update_matrix)
    signals = []
    outputs = []
    for i in eachindex(Nodes)
        if sum(update_matrix[:,i]) == 0
            push!(signals, i)
        end
        if sum(update_matrix[i,:]) == 0
            push!(outputs, i)
        end
    end
    signalNodes = Nodes[signals]
    outputNodes = Nodes[outputs]
    if repeat
        while length(signals) + length(outputs) > 0
            peripherals = vcat(signals, outputs)
            remaining = setdiff(1:length(Nodes), peripherals)
            #remove peripheral nodes from Nodes vector
            Nodes = Nodes[remaining]
            #remove peripheral rows and columns
            update_matrix = update_matrix[remaining, remaining]
            signals = []
            outputs = []
            for i in eachindex(Nodes)
                if sum(update_matrix[:,i]) == 0
                    push!(signals, i)
                end
                if sum(update_matrix[i,:]) == 0
                    push!(outputs, i)
                end
            end
            if length(signals) != 0
                signalNodes = vcat(signalNodes, Nodes[signals])
            end
            if length(outputs) != 0
                outputNodes = vcat(outputNodes, Nodes[outputs])
            end
        end
    end
    return signalNodes, outputNodes
end

function defaultWeightsFunction(noise::Float64)
    function weightsFunction(randVec::Array{Float64,1})
        #normal random variable with mean 0 and variance noise
        randN = randn(length(randVec))*noise
        rVec = Float64[]
        for i in eachindex(randVec)
            if (abs(randVec[i]) != 1)
                push!(rVec, min(max(randVec[i] + randN[i], 0), 1))
            else
                push!(rVec, randVec[i])
            end
        end
        return rVec
    end
    return weightsFunction
end

function contWeightUp(noise::Float64, update_matrix::Union{Array{Int,2}, Array{Float64,2}})
    function weightsFunction(randVec::Array{Float64, 1})
        
    end
    return weightsFunction
end

function adjust_indices(indices, deleted)
    sort!(deleted)
    adjusted = Int[]

    for idx in indices
        ct = count(x -> x < idx, deleted)
        push!(adjusted, idx - ct)
    end

    return adjusted
end



function load_bn_from_json(filename::String)
    data = JSON.parsefile(filename)

    F = [Vector{Int}(row) for row in data["F"]]
    I = [Vector{Int}(row) for row in data["I"]]
    I = [i .+ 1 for i in I]  # Adjust indices to be 1-based
    degree = Vector{Int}(data["degree"])
    variables = Vector{String}(data["variables"])
    constants = Vector{String}(data["constants"])

    return F, I, length(F), degree, variables, constants
end




function dec2binvec(n, Nbits)
    dig = reverse(digits(n, base=2, pad = Nbits))
    return dig
end

function pyCheck()
    modules = ["numpy", "random", "itertools", "json"]
    counter = 0
    for i in modules
        try
            pyimport(i)
        catch e
            println("Python module $i not found.")
            counter +=1
        end
    end
    if counter > 0
        println("Please install the missing modules.")
        return false
    end
    return true
end

function JSD(p::Vector{Float64},q::Vector{Float64}; b::Int = 2)
    mid = (p + q)./2
    return (kldivergence(p,mid, b) + kldivergence(q, mid, b))./2
end

function JSD(x::Matrix{Float64}; b::Int = 2)
    N = size(x, 1)
    jsd = Dict{Tuple{Int, Int}, Float64}()
    for i in 1:(N-1)
        for j in (i+1):N
            jsd[(i,j)] = JSD(x[i,:], x[j,:]; b = b)
        end
    end
    return jsd
end