# function text_to_BN(folder::String, textfile::String; 
#     separator_var_func::String="=", 
#     original_not::String="NOT", 
#     original_and::String="AND", 
#     original_or::String="OR", 
#     new_not::String="not", 
#     new_and::String="and", 
#     new_or::String="or", 
#     max_degree::Int=15, 
#     TREATMENT_OF_CONSTANTS::Int=1, 
#     max_N::Int=10000)
# """
# Translates a text file describing gene regulation into a Boolean network model.
# See original docstring for details.
# """

# # Read and preprocess the file
#     text = read(joinpath(folder, textfile), String)
#     text = replace(text, "\t" => " ", "(" => " ( ", ")" => " ) ")
#     tvec = split(text, "\n")
#     tvec = filter(!isempty, tvec)

#     n = length(tvec)
#     @assert n <= max_N "n=$n > max_N=$max_N"

#     # Extract variable names
#     var = Vector{String}(undef, n)
#     for i in 1:n
#     idx = findfirst(separator_var_func, tvec[i])
#     var[i] = replace(tvec[i][1:idx-1], " " => "")
#     end

#     # Identify constants and variables
#     constants_and_variables = String[]
#     for line in tvec
#     linesplit = split(line)
#     for el in linesplit
#     if el ∉ ["(", ")", "+", "*", "1", separator_var_func, original_not, original_and, original_or, "", " "]
#     push!(constants_and_variables, el)
#     end
#     end
#     end
#     constants = setdiff(unique(constants_and_variables), var)

#     # Create dictionary for variable and constant mapping
#     dict_variables_and_constants = Dict(
#     original_not => new_not,
#     original_and => new_and,
#     original_or => new_or
#     )
#     for (i, v) in enumerate(var)
#     dict_variables_and_constants[v] = "x[$i]"
#     end
#     for (i, c) in enumerate(constants)
#     dict_variables_and_constants[c] = "x[$(length(var) + i)]"
#     end

#     # Replace variables and operators in tvec
#     for (i, line) in enumerate(tvec)
#     linesplit = split(line)
#     for (ii, el) in enumerate(linesplit)
#     if el ∉ ["(", ")", "+", "*", "1", separator_var_func, new_not, new_and, new_or, "", " "]
#     linesplit[ii] = get(dict_variables_and_constants, el, el)
#     end
#     end
#     tvec[i] = join(linesplit, " ")
#     end

#     # Extract update rules
#     for i in 1:n
#     idx = findfirst(separator_var_func, tvec[i])
#     tvec[i] = tvec[i][idx + length(separator_var_func):end]
#     end

#     # Create I: list of regulator indices
#     I = Vector{Vector{Int}}()
#     tvec_mod = Vector{String}()
#     for i in 1:n
#     indices_open = findall("[", tvec[i])
#     indices_end = findall("]", tvec[i])
#     indices = sort(unique(parse.(Int, [tvec[i][begin+1:end-1] for (begin, end) in zip(indices_open, indices_end)])))
#     push!(I, indices)
#     dict_dummy = Dict(idx => j for (j, idx) in enumerate(indices))
#     tvec_dummy = tvec[i]
#     for el in sort(indices)
#     tvec_dummy = replace(tvec_dummy, "[$el]" => "[$(dict_dummy[el])]")
#     end
#     push!(tvec_mod, tvec_dummy)
#     end

#     degree = length.(I)

#     # Compute truth tables (F) using a custom parser
#     F = Vector{Vector{Int}}()
#     for i in 1:n
#     f = Int[]
#     if degree[i] <= max_degree
#     X = collect(Iterators.product(fill(0:1, degree[i])...))
#     expr = tvec_mod[i]
#     for x in X
#     result = evaluate_boolean(expr, collect(x))
#     push!(f, result)
#     end
#     end
#     push!(F, f)
#     end

#     # Handle constants
#     if TREATMENT_OF_CONSTANTS == 1
#     for i in 1:length(constants)
#     push!(F, [0, 1])
#     push!(I, [length(var) + i])
#     push!(degree, 1)
#     end
#     end
#     @assert TREATMENT_OF_CONSTANTS in [0, 1] "TREATMENT_OF_CONSTANTS must be 0 or 1. Option 2 is not implemented."

#     return F, I, degree, var, constants
# end

# # Helper function to find all indices of a substring
# function findall(substr::String, str::String)
# indices = Int[]
# start = 1
# while true
# idx = findnext(substr, str, start)
# if idx === nothing
# break
# end
# push!(indices, first(idx))
# start = last(idx) + 1
# end
# return indices
# end

# # Custom Boolean expression parser and evaluator
# function evaluate_boolean(expr::String, x::Vector{Int})
# # Tokenize the expression
# tokens = split(replace(expr, "(" => " ( ", ")" => " ) "), r"\s+"; keepempty=false)

# # Function to evaluate tokens recursively
# function eval_tokens(tokens::Vector{String}, pos::Int=1, x::Vector{Int}=x)
# result = nothing
# op = nothing
# i = pos

# while i <= length(tokens)
# token = tokens[i]

# if token == "("
# # Start a new subexpression
# sub_result, new_i = eval_tokens(tokens, i + 1, x)
# i = new_i
# if result === nothing
#     result = sub_result
# elseif op == "not"
#     result = 1 - result
#     op = nothing
# elseif op == "and"
#     result = result & sub_result
#     op = nothing
# elseif op == "or"
#     result = result | sub_result
#     op = nothing
# end

# elseif token == ")"
# # End of subexpression
# return result, i + 1

# elseif token == "not"
# op = "not"
# i += 1

# elseif token == "and"
# op = "and"
# i += 1

# elseif token == "or"
# op = "or"
# i += 1

# elseif occursin(r"^\[(\d+)\]$", token)
# # Variable reference, e.g., [1]
# idx = parse(Int, match(r"\[(\d+)\]", token).captures[1])
# val = x[idx]

# if result === nothing
#     result = val
# elseif op == "not"
#     result = 1 - val
#     op = nothing
# elseif op == "and"
#     result = result & val
#     op = nothing
# elseif op == "or"
#     result = result | val
#     op = nothing
# end
# i += 1

# else
# error("Invalid token: $token")
# end
# end

# return result, i
# end

# # Evaluate the expression
# result, _ = eval_tokens(tokens, 1, x)
# return result % 2  # Ensure binary output
# end