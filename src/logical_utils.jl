### This file has accompany functions for logical rule processing and simulation. 
### List of functions here:
# read_boolean_rules -> dataframe with target and rules
# tokenize_expression -> get token out of a rule, which basically means the node names, brackets and logical operations. Processes the node names by default and can get a map of unprocessed version too. - Get a dict instead?
# infix_to_postfix -> convert the tokens to postfix
# postfix_to_dnf -> convert postfix to a list of list such that each sublist is exclusively an and function (can we come up with a way to do this for ising rule too?)
# dnf_to_truthtable -> evaluate the dnf to get the truth table pretty straightforward
# f_from_expr -> For a given rule/expr, tokenize it, get postfix, convert to dnf and extract truthtable. Returns truthtable, dnf and vars.
# f_from_file -> get the list of truth tables from a file showing Boolean rules
# parse_boolean_network -> Given a boolean network, returns F (list of truth tables), I, N, degrees, variables, constants
# create_dnf_evaluator -> for a dnf and a node_to_ids dictionary, generate a function that takes in a state and evaluates the output for that dnf - I assume it to be faster than truth table based evaluation, because truthtable lookups require multiple bin2dec and dec2bin calculations, while this won't. It should also be scalable, since the truthtable method encodes states as integers and larger network states turn into huge decimals. Not sure if that is true though.
# getNodeFunctions -> from a boolean network file, get update functions for each node. guess I could do it from the truth tables too, but doesn't matter either way. Once a network, so overhead is very small. Would be even smaller if I have a class.
# is_monotonic -> given a truth table, evaluates if the interaction from each source to the target is monotonic or not
# boolean_to_topo -> Takes a boolean network file, calculates the turth tables and evaluates the nature of the interactions, stores them as a topo file
# F_to_topo -> Gets a topo file out of a list of truth tables
###


# -------------------------------
# 1. Read Boolean Rules
# -------------------------------
function read_boolean_rules(filename::String; sep::String=" = ")
    lines = readlines(filename)
    # remove comments (#) and blank lines
    rules_raw = filter(l -> !(startswith(l, "#") || isempty(strip(l))), lines)

    targets = String[]
    rules = String[]

    for line in rules_raw
        parts = split(line, sep, limit=2)
        if length(parts) == 2
            target = strip(parts[1])
            rule = strip(parts[2])
            push!(targets, target)
            push!(rules, rule)
        end
    end

    return DataFrame(target = targets, rule = rules)
end

# -------------------------------
# 2. Tokenize Expression
# -------------------------------
function tokenize_expression(expr::String, getOrig::Bool)
    expr = strip(expr)
    if (!getOrig)
        expr = uppercase(expr)
    end
    expr = replace(expr, r"\s+" => " ")
    expr = replace(expr, "(" => " ( ")
    expr = replace(expr, ")" => " ) ")
    tokens = split(expr, r"\s+")
    clean_tokens = String[]

    for t in tokens
        if t in ["(", ")"]
            push!(clean_tokens, t)
        else
            if (!getOrig)
                t_clean = replace(t, r"\W+" => "")
                t_clean = replace(t_clean, "_" => "")
            else
                t_clean = t
            end
            if !isempty(t_clean)
                push!(clean_tokens, t_clean)
            end
        end
    end

    return clean_tokens
end

# -------------------------------
# 3. Infix → Postfix (Shunting Yard)
# -------------------------------
function infix_to_postfix(tokens::Vector{String})
    precedence = Dict("NOT" => 3, "AND" => 2, "OR" => 1, "(" => 0)
    operators = Set(["NOT", "AND", "OR"])
    output = String[]
    stack = String[]

    for token in tokens
        if token in operators
            while !isempty(stack) && (last(stack) in operators) &&
                  precedence[last(stack)] ≥ precedence[token]
                push!(output, pop!(stack))
            end
            push!(stack, token)
        elseif token == "("
            push!(stack, token)
        elseif token == ")"
            while !isempty(stack) && last(stack) != "("
                push!(output, pop!(stack))
            end
            if !isempty(stack) && last(stack) == "("
                pop!(stack)
            end
        else
            push!(output, token)
        end
    end

    while !isempty(stack)
        push!(output, pop!(stack))
    end

    return output
end

# -------------------------------
# 4. Postfix → DNF (list of list of literals)
# -------------------------------
function postfix_to_dnf(postfix::Vector{String})
    stack = Vector{Any}()

    for token in postfix
        if token == "NOT"
            if !isempty(stack)
                operand = pop!(stack)
                # operand is a DNF: [[A, B], [C]] means (A AND B) OR C
                # NOT((A AND B) OR C) = NOT(A AND B) AND NOT(C)
                #                     = (NOT A OR NOT B) AND NOT C
                # Convert back to DNF: (NOT A AND NOT C) OR (NOT B AND NOT C)
                
                # Step 1: Convert each DNF term to a CNF clause by negating literals
                # Term [A, B] (meaning A AND B) becomes [~A, ~B] (meaning ~A OR ~B)
                cnf_clauses = []
                for term in operand
                    clause = [startswith(lit, "~") ? lit[2:end] : "~" * lit for lit in term]
                    push!(cnf_clauses, clause)
                end
                
                # Step 2: Convert CNF to DNF by distributing ANDs over ORs
                # We need the cross product of all clauses
                if isempty(cnf_clauses)
                    result = [[]]  # Empty DNF = True
                else
                    # Start with first clause - each literal becomes a term
                    result = [[lit] for lit in cnf_clauses[1]]
                    
                    # For each additional clause, combine with existing terms
                    for clause in cnf_clauses[2:end]
                        new_result = []
                        for existing_term in result
                            for lit in clause
                                # Combine existing_term AND lit
                                push!(new_result, vcat(existing_term, [lit]))
                            end
                        end
                        result = new_result
                    end
                end
                
                push!(stack, result)
            end
        elseif token == "AND"
            if length(stack) ≥ 2
                op2 = pop!(stack)
                op1 = pop!(stack)
                result = Any[]
                for term1 in op1, term2 in op2
                    push!(result, vcat(term1, term2))
                end
                push!(stack, result)
            end
        elseif token == "OR"
            if length(stack) ≥ 2
                op2 = pop!(stack)
                op1 = pop!(stack)
                push!(stack, vcat(op1, op2))
            end
        else
            push!(stack, [[token]])
        end
    end

    return isempty(stack) ? Vector{String}[] : [Vector{String}(x) for x in stack[end]]
end

# -------------------------------
# 5. DNF → Truth Table
# -------------------------------
"""
    dnf_to_truthtable(dnf_terms::Vector{Vector{String}}) -> DataFrame

Given a DNF (list of terms, each a list of literals),
returns the full truth table with an Output column.
"""
function dnf_to_truthtable(dnf_terms::Vector{Vector{String}})
    # collect unique variables
    vars = Set{String}()
    for term in dnf_terms, lit in term
        varname = startswith(lit, "~") ? lit[2:end] : lit
        push!(vars, varname)
    end
    varlist = sort(collect(vars))

    nvars = length(varlist)
    nrows = 2^nvars
    df = DataFrame([var => falses(nrows) for var in varlist])

    # fill truth assignments
    for (i, row) in enumerate(eachrow(df))
        for (j, var) in enumerate(varlist)
            row[var] = (i-1) >> (nvars - j) & 1 == 1
        end
    end

    output = falses(nrows)

    # Evaluate DNF for each row
    for (i, row) in enumerate(eachrow(df))
        # OR over terms
        term_vals = Bool[]
        for term in dnf_terms
            lit_vals = Bool[]
            for lit in term
                if startswith(lit, "~")
                    var = lit[2:end]
                    push!(lit_vals, !row[var])
                else
                    push!(lit_vals, row[lit])
                end
            end
            push!(term_vals, all(lit_vals))
        end
        output[i] = any(term_vals)
    end

    df[!, "Output"] = output
    return df
end

# Example usage:
# tokens = tokenize_expression("A AND (NOT B OR C)")
# postfix = infix_to_postfix(tokens)
# dnf = postfix_to_dnf(postfix)
# tt = dnf_to_truthtable(dnf)
# println(tt)


"""
    f_from_expression(expr::String) -> Tuple{Vector{Int}, Vector{String}}

Parse a Boolean expression and return (truth_table, variable_names).
Uses the tokenize/postfix/DNF pipeline.
"""
function f_from_expression(expr::String)
    tokens = tokenize_expression(expr, false)
    tokenOrig = tokenize_expression(expr, true)
    postfix = infix_to_postfix(tokens)
    dnf = postfix_to_dnf(postfix)
    tt_df = dnf_to_truthtable(dnf)
    tokens = vcat([tokens, ["~"*x for x in tokens]]...)
    tokenOrig = vcat([tokenOrig, ["~"*x for x in tokenOrig]]...)
    tokenDict = Dict(tokens .=> tokenOrig)
    # Extract variable names (all columns except "Output")
    vars = [name for name in names(tt_df) if name != "Output"]
    vars = [tokenDict[i] for i in vars]
    # Extract truth table values
    truth_table = Int.(tt_df.Output)
    dnf = [[tokenDict[i] for i in x] for x in dnf]
    return truth_table, vars, dnf
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

    table = Dict{String, Tuple{Vector{Int}, Vector{String}, Vector{Vector{String}}}}()
    for line in fl
        name_expr = split(line, "=")
        name = string(strip(name_expr[1]))
        expr = string(name_expr[2])
        f, vars, dnf = f_from_expression(expr)
        table[name] = (f, vars, dnf)
    end
    return table
end



"""
    parse_boolean_network(rules_file::String) -> Tuple

Complete parsing pipeline: reads Boolean rules and returns simulation-ready format.

Arguments:
- `rules_file`: Path to file with Boolean rules (format: "NodeName = Boolean Expression")

Returns:
- `F`: Vector of truth tables for each node
- `I`: Vector of regulator index vectors
- `N`: Number of nodes
- `node_names`: Vector of node names
- `node_to_idx`: Dictionary mapping node names to indices
"""
function parse_boolean_network(rules_file::String)
    # Use existing parser
    truth_table_dict = f_from_file(rules_file)
    
    # Create node name to index mapping
    node_names = sort(unique(vcat([x[2] for x in collect(values(truth_table_dict))]...)))
    nd2 = collect(keys(truth_table_dict))
    node_names = sort(unique(vcat(node_names, nd2)))
    N = length(node_names)
    node_to_idx = Dict(name => i for (i, name) in enumerate(node_names))
    
    # Initialize F and I
    F = Vector{Vector{Int}}(undef, N)
    I = Vector{Vector{Int}}(undef, N)
    variables = Vector{String}()
    constants = Vector{String}()
    # Fill F and I for each node
    for (i, node_name) in enumerate(node_names)
        truth_table, regulators = get(truth_table_dict, node_name, (nothing, nothing))
        
        # Store truth table
        
        
        # Convert regulator names to indices
        if !isnothing(regulators)
            F[i] = truth_table
            regulator_indices = [node_to_idx[reg] for reg in regulators]
            I[i] = regulator_indices
            if (I[i] == [i] && F[i] == [0,1])
                push!(constants, node_name)
            else
                push!(variables, node_name)
            end
        else
            F[i] = [0,1]
            I[i] = [i]
            push!(constants, node_name)
        end
    end
    degrees = [length(x) for x in I]
    
    
    return F, I, N, node_names, degrees, variables, constants
end


"""
    create_dnf_evaluator(dnf::Vector{Vector{String}}, 
                        node_to_idx::Dict{String,Int},
                        token_map::Dict{String,String}) -> Function

Create a function that evaluates a DNF expression.

DNF structure: OR of terms, where each term is AND of literals.
- Literal "A" means variable A must be 1
- Literal "~A" means variable A must be 0

Arguments:
- `dnf`: DNF as list of terms, each term is list of literals
- `node_to_idx`: Mapping from node names to state vector indices
- `token_map`: Mapping from uppercase tokens to original case tokens

Returns:
- Function that takes state vector and returns 0 or 1
"""
function create_dnf_evaluator(dnf::Vector{Vector{String}}, 
                              node_to_idx::Dict{String,Int})
    
    # Pre-process DNF to avoid lookups during evaluation
    # Convert to: Vector of (is_negated, node_index) tuples
    processed_dnf = Vector{Vector{Tuple{Bool,Int}}}()
    
    for term in dnf
        processed_term = Tuple{Bool,Int}[]
        for literal in term
            is_negated = startswith(literal, "~")
            varname = is_negated ? literal[2:end] : literal
            
            # Get original token name
            original_name = varname
            
            # Get node index
            if haskey(node_to_idx, original_name)
                node_idx = node_to_idx[original_name]
                push!(processed_term, (is_negated, node_idx))
            else
                @warn "Variable $original_name not found in node mapping"
            end
        end
        push!(processed_dnf, processed_term)
    end
    
    # Return evaluation function
    return function evaluate_dnf(x::Vector{Int})
        # OR over terms
        for term in processed_dnf
            # AND over literals in term
            term_satisfied = true
            for (is_negated, node_idx) in term
                node_value = x[node_idx]
                if is_negated
                    if node_value == 1
                        term_satisfied = false
                        break
                    end
                else
                    if node_value == 0
                        term_satisfied = false
                        break
                    end
                end
            end
            
            if term_satisfied
                return 1
            end
        end
        
        return 0
    end
end

function getNodeFunctions(rules_file::String)
    truth_table_dict = f_from_file(rules_file)
    
    # Create node name to index mapping
    node_names = sort(unique(vcat([x[2] for x in collect(values(truth_table_dict))]...)))
    nd2 = collect(keys(truth_table_dict))
    node_names = sort(unique(vcat(node_names, nd2)))
    N = length(node_names)
    node_to_idx = Dict(name => i for (i, name) in enumerate(node_names))
    
    # Initialize F and I
    F = Vector{Function}(undef, N)
    I = Vector{Vector{Int}}(undef, N)
    variables = Vector{String}()
    constants = Vector{String}()
    # Fill F and I for each node
    for (i, node_name) in enumerate(node_names)
        truth_table, regulators, dnf = get(truth_table_dict, node_name, (nothing, nothing, nothing))
        
        # Convert regulator names to indices
        if !isnothing(regulators)
            F[i] = create_dnf_evaluator(dnf, node_to_idx)
            regulator_indices = [node_to_idx[reg] for reg in regulators]
            I[i] = regulator_indices
            if (I[i] == [i] && truth_table == [0,1])
                push!(constants, node_name)
            else
                push!(variables, node_name)
            end
        else
            F[i] = create_dnf_evaluator([[node_name]], node_to_idx)
            I[i] = [i]
            push!(constants, node_name)
        end
    end
    degrees = [length(x) for x in I]

    return F, I, N, node_names, degrees, variables, constants

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

