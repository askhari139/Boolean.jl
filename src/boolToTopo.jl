
"""
    f_from_expression(expr::String) -> Tuple{Vector{Int}, Vector{String}}

Convert a Boolean expression into its truth table representation (vector of 0s and 1s) 
and a list of input variables.
"""
function f_from_expression(expr::String)
    expr = strip(expr)
    expr = replace(expr, "AND" => "&", "OR" => "|", "NOT" => "!")  # Julia uses `!` for NOT
    expr = replace(expr, "(" => " ( ", ")" => " ) ")
    expr_split = split(expr, " ")
    
    var = String[]
    dict_var = Dict{String, String}()
    n_var = 0
    id_var = 1

    for (i, el) in enumerate(expr_split)
        if !(el in ["", " ", "(", ")", "and", "or", "not", "&", "|", "!", "+", "-", "*", "%", ">", ">=", "==", "<=", "<"]) && !occursin(r"^\d+$", el)
            if !haskey(dict_var, el)
                dict_var[el] = "x[$id_var]"
                push!(var, el)
                n_var += 1
                id_var += 1
            end
            expr_split[i] = dict_var[el]
        end
    end

    expr_reconstructed = join(expr_split, " ")
    println(expr_reconstructed)
    F = Int[]
    for x in IterTools.product(([0,1] for _ in 1:n_var)...)
        x = reverse(collect(x))
        push!(F, Int(eval(Meta.parse(expr_reconstructed))))
    end

    return F, var
end

"""
    is_monotonic(F::Vector{Int}; GET_DETAILS::Bool=false) -> Bool or Tuple

Check whether a Boolean function `F` is monotonic.
If `GET_DETAILS=true`, returns a tuple (is_monotonic, regulation_types).
"""
function is_monotonic(F::Vector{Int}; GET_DETAILS::Bool=false)
    n = Int(log2(length(F)))
    F = collect(F)
    monotonic = String[]
    for i in 1:n
        dummy_add = 2^(n - i)
        dummy = [(j % 2^(n - i + 1)) ÷ dummy_add for j in 0:(2^n - 1)]
        diff = [F[j+1] for j in 0:(2^n - 1) if dummy[j+1] == 1] .- [F[j+1] for j in 0:(2^n - 1) if dummy[j+1] == 0]
        min_diff = minimum(diff)
        max_diff = maximum(diff)
        if min_diff == 0 && max_diff == 0
            push!(monotonic, "0")
        elseif min_diff == -1 && max_diff == 1
            push!(monotonic, "0")
        elseif min_diff ≥ 0 && max_diff == 1
            push!(monotonic, "1")
        elseif min_diff == -1 && max_diff ≤ 0
            push!(monotonic, "2")
        else
            push!(monotonic, "0")
        end
    end
    return GET_DETAILS ? (!("0" in monotonic), monotonic) : !("0" in monotonic)
end

"""
    f_from_file(filename::String) -> Dict{String, Tuple{Vector{Int}, Vector{String}}}

Reads Boolean rules from a file and returns a mapping from node name to (truth table, variable list).
"""
function f_from_file(filename::String)
    fl = String[]
    open(filename, "r") do f
        for line in eachline(f)
            if !startswith(line, "#")
                push!(fl, strip(line))
            end
        end
    end

    table = Dict{String, Tuple{Vector{Int}, Vector{String}}}()
    for line in fl
        name_expr = split(line, "=")
        name = string(strip(name_expr[1]))
        expr = string(name_expr[2])
        f, vars = f_from_expression(expr)
        table[name] = (f, vars)
    end
    return table
end

"""
    boolean_to_topo(input_file::String, output_file::String)

Generates a .topo file from a Boolean rules file describing regulatory types (1 = activation, 2 = inhibition).
"""
function boolean_to_topo(input_file::String)
    TruthTable = f_from_file(input_file)
    topoFile = ["Source Target Type"]
    # remove .txt from input file name 
    output_file = replace(input_file, r"\.txt$" => ".topo")
    for (key, (f, vars)) in TruthTable
        _, regulation_types = is_monotonic(f, GET_DETAILS=true)
        for (i, node) in enumerate(vars)
            if regulation_types[i] != "0"
                push!(topoFile, "$node $key $(regulation_types[i])")
            end
        end
    end

    open(output_file, "w") do f
        for line in topoFile
            println(f, line)
        end
    end
end

"""
    F_to_topo(TruthTable::Dict{String, Tuple{Vector{Int}, Vector{String}}}, output_file::String)

Write topology to file given a preloaded truth table (instead of a file).
"""
function F_to_topo(TruthTable::Dict{String, Tuple{Vector{Int}, Vector{String}}}, output_file::String)
    topoFile = ["Source Target Type"]

    for (key, (f, vars)) in TruthTable
        _, regulation_types = is_monotonic(f, GET_DETAILS=true)
        for (i, node) in enumerate(vars)
            if regulation_types[i] != "0"
                push!(topoFile, "$node $key $(regulation_types[i])")
            end
        end
    end

    open(output_file * ".topo", "w") do f
        for line in topoFile
            println(f, line)
        end
    end
end
