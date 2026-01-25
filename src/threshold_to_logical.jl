#### The output from these boolean rules still does not match that of the ising simulations. 
function nchoosek(n,k)
    factorial(n)/(factorial(k)*factorial(n-k))
end

"""
    calculate_threshold_conditions_general(
        P::Int,
        N::Int, 
        SR::Real,
        theta::Real;
        w_pos::Real = 1.0,
        w_neg::Real = -1.0,
        state_values::Tuple = (0, 1)
    ) -> Vector{NamedTuple}

Calculate threshold-based update conditions for a node with P positive regulators,
N negative regulators, and self-regulation weight SR.

The update rule is:
    s_i(t+1) = s_high if Σ J_ji * s_j > theta
    s_i(t+1) = s_low  if Σ J_ji * s_j < theta
    s_i(t+1) = s_i(t) if Σ J_ji * s_j == theta

# Arguments
- `P::Int`: Number of positive (activating) regulators
- `N::Int`: Number of negative (inhibiting) regulators
- `SR::Real`: Self-regulation weight (default: 0 for no self-regulation; 1 for self activation and -1 for self inhibition)
- `theta::Real`: Threshold value
- `w_pos::Real`: Weight for positive regulators (default: 1.0)
- `w_neg::Real`: Weight for negative regulators (default: -1.0)
- `state_values::Tuple`: (s_low, s_high) possible state values (default: (0, 1))

# Returns
Vector of NamedTuples with fields:
- `P_plus_N_minus`: Value of P+ + N- (active activators + inactive inhibitors)
- `sum_value`: The actual weighted sum (excluding self-regulation)
- `si_low_output`: Output when current state is s_low
- `si_high_output`: Output when current state is s_high
- `si_condition`: Symbol indicating if output depends on current state
  - `:none` - output same regardless of s_i
  - `:positive` - output depends on s_i being high
  - `:negative` - output depends on s_i being low
  - `:required` - output equals s_i (at threshold)
- `output`: The output value (if si_condition == :none)

# Example
```julia
# Standard Ising model with 3 activators, 2 inhibitors, self-activation
conditions = calculate_threshold_conditions_general(3, 2, 1.0, 2.5)

# Boolean network with (0,1) states and theta = 0.5*(P+N)
conditions = calculate_threshold_conditions_general(
    3, 2, 0.0, 2.5, 
    state_values=(0, 1)
)
```
"""
function calculate_threshold_conditions_general(
    P::Int,
    N::Int, 
    SR::Real;
    theta::Real = (P-N)/2,
    w_pos::Real = 1.0,
    w_neg::Real = -1.0,
    state_values::Tuple = (0, 1)
)
    println(theta)
    s_low, s_high = state_values
    
    # Store all conditions
    conditions = NamedTuple[]
    
    # Iterate through all possible combinations of active regulators
    for P_plus in 0:P  # number of active activators
        for N_plus in 0:N  # number of inactive inhibitors
            
            # Calculate external contribution (excluding self-regulation)
            sum_ext = P_plus * w_pos * s_high +           # active activators
                     (P - P_plus) * w_pos * s_low +       # inactive activators
                     N_plus * w_neg * s_high +     # active inhibitors
                     (N - N_plus) * w_neg * s_low              # inactive inhibitors
            
            # Determine output for each possible current state
            sum_with_si_low = sum_ext + SR * s_low
            sum_with_si_high = sum_ext + SR * s_high
            
            # Determine outputs
            if sum_with_si_low > theta
                output_si_low = s_high
            elseif sum_with_si_low < theta
                output_si_low = s_low
            else
                output_si_low = s_low  # stays at current state
            end
            
            if sum_with_si_high > theta
                output_si_high = s_high
            elseif sum_with_si_high < theta
                output_si_high = s_low
            else
                output_si_high = s_high  # stays at current state
            end
            
            # Determine si_condition and output
            si_condition = :none
            output = nothing
            
            if sum_with_si_low == theta && sum_with_si_high == theta
                # Both at threshold - output always equals input
                si_condition = :positive
                output = nothing
            elseif output_si_low == output_si_high
                # Output doesn't depend on current state
                si_condition = :none
                output = output_si_low
            elseif sum_with_si_low == theta
                # Only low state at threshold
                si_condition = :negative
                output = output_si_high  # definite output when si = s_high
            elseif sum_with_si_high == theta
                # Only high state at threshold
                si_condition = :positive
                output = output_si_low  # definite output when si = s_low
            elseif output_si_low == s_low && output_si_high == s_high
                # Output follows current state (sum straddles threshold)
                si_condition = :positive  # need si=high to get high output
                output = nothing
            elseif output_si_low == s_high && output_si_high == s_low
                # Output opposite of current state (rare case)
                si_condition = :negative  # need si=low to get high output
                output = nothing
            end
            
            push!(conditions, (
                P_plus_N_plus = P_plus + N_plus,
                P_plus = P_plus,
                N_plus = N_plus,
                n_terms = nchoosek(P, P_plus)*nchoosek(N, N_plus),
                sum_value = sum_ext,
                si_low_output = output_si_low,
                si_high_output = output_si_high,
                si_condition = si_condition,
                output = output
            ))
        end
    end
    
    return conditions
end

# function test()
#     println("Revise working")
# end
"""
    parse_literal(literal::String) -> Tuple{String, Bool}

Parse a literal into (variable_name, is_positive).
Returns ("A", true) for "A" and ("A", false) for "~A".
"""
function parse_literal(literal::String)
    if startswith(literal, "~")
        return (literal[2:end], false)
    else
        return (literal, true)
    end
end


"""
    clause_to_dict(clause::Vector{String}) -> Dict{String, Bool}

Convert a clause to a dictionary mapping variable -> polarity.
"""
function clause_to_dict(clause::Vector{String})
    result = Dict{String, Bool}()
    for lit in clause
        var, pol = parse_literal(lit)
        result[var] = pol
    end
    return result
end


"""
    dict_to_clause(d::Dict{String, Bool}) -> Vector{String}

Convert a variable->polarity dictionary back to a clause.
"""
function dict_to_clause(d::Dict{String, Bool})
    clause = String[]
    for (var, pol) in sort(collect(d), by=x->x[1])  # Sort for consistency
        if pol
            push!(clause, var)
        else
            push!(clause, "~" * var)
        end
    end
    return clause
end


"""
    subsumes(clause1::Vector{String}, clause2::Vector{String}) -> Bool

Check if clause1 subsumes clause2 (i.e., clause1 ⊆ clause2).
If clause1 ⊆ clause2, then clause2 is redundant: (A) OR (A AND B) = A
"""
function subsumes(clause1::Vector{String}, clause2::Vector{String})
    set1 = Set(clause1)
    set2 = Set(clause2)
    return issubset(set1, set2)
end


"""
    remove_subsumed_clauses(dnf::Vector{Vector{String}}) -> Vector{Vector{String}}

Remove clauses that are subsumed by other clauses.
"""
function remove_subsumed_clauses(dnf::Vector{Vector{String}})
    n = length(dnf)
    keep = trues(n)
    
    for i in 1:n
        if !keep[i]
            continue
        end
        for j in 1:n
            if i == j || !keep[j]
                continue
            end
            # If clause i subsumes clause j, remove j
            if subsumes(dnf[i], dnf[j])
                keep[j] = false
            end
        end
    end
    
    return [dnf[i] for i in 1:n if keep[i]]
end


"""
    find_common_and_differing(clauses::Vector{Vector{String}}) 
        -> Tuple{Dict{String,Bool}, Set{String}}

Given a set of clauses, find:
1. Common literals (variables with same polarity in all clauses)
2. Differing variables (variables that have different polarities across clauses)

Returns (common_dict, differing_vars)
"""
function find_common_and_differing(clauses::Vector{Vector{String}})
    if isempty(clauses)
        return (Dict{String,Bool}(), Set{String}())
    end
    
    # Convert all clauses to dicts
    clause_dicts = [clause_to_dict(c) for c in clauses]
    
    # Find all variables that appear
    all_vars = Set{String}()
    for d in clause_dicts
        union!(all_vars, keys(d))
    end
    
    # Check each variable
    common = Dict{String, Bool}()
    differing = Set{String}()
    
    for var in all_vars
        # Check if this variable appears in all clauses
        appears_in_all = all(haskey(d, var) for d in clause_dicts)
        
        if !appears_in_all
            # Variable doesn't appear in all clauses, can't simplify
            continue
        end
        
        # Check if polarity is same in all clauses
        polarities = [d[var] for d in clause_dicts]
        if all(p == polarities[1] for p in polarities)
            # Same polarity in all - this is a common literal
            common[var] = polarities[1]
        else
            # Different polarities - this is a differing variable
            push!(differing, var)
        end
    end
    
    return (common, differing)
end


"""
    check_complete_coverage(clauses::Vector{Vector{String}}, 
                           differing_vars::Set{String}) -> Bool

Check if the clauses provide complete Boolean coverage over the differing variables.
For k differing variables, we need 2^k clauses covering all combinations.
"""
function check_complete_coverage(clauses::Vector{Vector{String}}, 
                                differing_vars::Set{String})
    if isempty(differing_vars)
        return false  # No variables to absorb
    end
    
    k = length(differing_vars)
    expected_count = 2^k
    
    if length(clauses) != expected_count
        return false
    end
    
    # Extract the polarities of differing variables from each clause
    differing_list = sort(collect(differing_vars))
    patterns = Set{Vector{Bool}}()
    
    for clause in clauses
        d = clause_to_dict(clause)
        pattern = Bool[]
        
        # Check if all differing vars are present
        for var in differing_list
            if !haskey(d, var)
                return false  # Variable missing
            end
            push!(pattern, d[var])
        end
        
        push!(patterns, pattern)
    end
    
    # Check if we have all 2^k unique patterns
    return length(patterns) == expected_count
end


"""
    try_absorb_group(clauses::Vector{Vector{String}}) 
        -> Union{Vector{String}, Nothing}

Try to absorb a group of clauses into a simpler clause.
Returns the absorbed clause if successful, nothing otherwise.
"""
function try_absorb_group(clauses::Vector{Vector{String}})
    if length(clauses) < 2
        return nothing
    end
    
    common, differing = find_common_and_differing(clauses)
    
    if isempty(differing)
        # All clauses are identical (shouldn't happen after deduplication)
        return nothing
    end
    
    # Check if we have complete coverage
    if check_complete_coverage(clauses, differing)
        # We can absorb! Return just the common literals
        return dict_to_clause(common)
    end
    
    return nothing
end

"""
    simplify_dnf(dnf::Vector{Vector{String}}; max_iterations::Int=10) 
        -> Vector{Vector{String}}

Simplify a DNF expression by applying absorption and subsumption laws.

The function iteratively:
1. Removes duplicate clauses
2. Removes subsumed clauses (if A ⊆ B, remove B)
3. Applies absorption law: groups clauses with common literals and checks
   if they cover all combinations of differing variables

# Arguments
- `dnf`: DNF as vector of clauses (each clause is vector of literals)
- `max_iterations`: Maximum number of simplification passes

# Returns
Simplified DNF

# Example
```julia
# (A AND ~B AND ~C) OR (A AND B AND ~C) → (A AND ~C)
dnf = [["A", "~B", "~C"], ["A", "B", "~C"]]
simplified = simplify_dnf(dnf)  # Returns [["A", "~C"]]

# (A AND ~B AND ~C) OR (A AND B AND ~C) OR (A AND ~B AND C) OR (A AND B AND C) → [A]
dnf = [["A", "~B", "~C"], ["A", "B", "~C"], ["A", "~B", "C"], ["A", "B", "C"]]
simplified = simplify_dnf(dnf)  # Returns [["A"]]
```
"""
function simplify_dnf(dnf::Vector{Vector{String}}; max_iterations::Int=10)
    if isempty(dnf)
        return dnf
    end
    
    current_dnf = deepcopy(dnf)
    
    for iteration in 1:max_iterations
        # Remove duplicates
        current_dnf = unique(current_dnf)
        
        # Remove subsumed clauses
        current_dnf = remove_subsumed_clauses(current_dnf)
        
        # Try to find groups that can be absorbed
        changed = false
        
        # Group clauses by number of literals (optimization)
        # Clauses that will be absorbed must have the same number of literals
        by_length = Dict{Int, Vector{Int}}()
        for (i, clause) in enumerate(current_dnf)
            len = length(clause)
            if !haskey(by_length, len)
                by_length[len] = Int[]
            end
            push!(by_length[len], i)
        end
        
        new_dnf = Vector{String}[]
        used = falses(length(current_dnf))
        
        # Try to find absorbable groups
        for (len, indices) in by_length
            if len == 0 || length(indices) < 2
                continue
            end
            
            # Try all possible subsets of size 2^k for k = 1, 2, 3, ...
            for k in 1:min(5, len)  # Limit k to avoid exponential blowup
                subset_size = 2^k
                if subset_size > length(indices)
                    break
                end
                
                # Try all combinations of subset_size clauses
                for combo in combinations(indices, subset_size)
                    if any(used[i] for i in combo)
                        continue
                    end
                    
                    clauses_subset = [current_dnf[i] for i in combo]
                    absorbed = try_absorb_group(clauses_subset)
                    
                    if absorbed !== nothing
                        # Successfully absorbed!
                        push!(new_dnf, absorbed)
                        for i in combo
                            used[i] = true
                        end
                        changed = true
                    end
                end
            end
        end
        
        # Add remaining unused clauses
        for i in eachindex(current_dnf)
            if !used[i]
                push!(new_dnf, current_dnf[i])
            end
        end
        
        current_dnf = new_dnf
        
        # If no changes, we're done
        if !changed
            break
        end
    end
    
    # Final cleanup
    current_dnf = unique(current_dnf)
    current_dnf = remove_subsumed_clauses(current_dnf)
    
    # Sort clauses for consistent output
    sort!(current_dnf, by = c -> (length(c), join(sort(c))))
    
    return current_dnf
end



"""
    interaction_matrix_to_dnf(
        interaction_matrix::Matrix,
        node_names::Vector{String},
        node_index::Int;
        theta::Union{Real, Function, Nothing} = nothing,
        state_values::Tuple = (0, 1),
        theta_mode::Symbol = :half_total
    ) -> Tuple{Vector{Vector{String}}, Vector{String}}

Convert interaction matrix to DNF (Disjunctive Normal Form) for a specific node.

# Arguments
- `interaction_matrix::Matrix`: N×N interaction matrix (positive = activation, negative = inhibition)
- `node_names::Vector{String}`: Names of all nodes
- `node_index::Int`: Index of target node to extract DNF for
- `theta::Union{Real, Function, Nothing}`: Threshold value
  - If `Real`: use this fixed threshold
  - If `Function`: called as theta(P, N, SR) to compute threshold
  - If `nothing`: determined by theta_mode
- `state_values::Tuple`: (s_low, s_high) possible state values (default: (0, 1))
- `theta_mode::Symbol`: How to calculate theta if not provided
  - `:half_total` - theta = (P + N) / 2 (standard Ising)
  - `:zero` - theta = 0
  - `:custom` - must provide theta argument

# Returns
- `dnf::Vector{Vector{String}}`: DNF representation (OR of AND clauses)
- `regulators::Vector{String}`: Names of regulators (for reference)

# Example
```julia
# 3-node network: A activates B, C inhibits B, B self-activates
J = [0 1 0; 0 1 -1; 0 0 0]
names = ["A", "B", "C"]
dnf, regs = interaction_matrix_to_dnf(J, names, 2)  # Get DNF for node B
```
"""
function interaction_matrix_to_dnf(
    interaction_matrix::Matrix,
    node_names::Vector{String},
    node_index::Int;
    theta::Union{Real, Function, Nothing} = nothing,
    state_values::Tuple = (-1, 1),
    theta_mode::Symbol = :zero
)
    s_low, s_high = state_values
    n_nodes = size(interaction_matrix, 1)
    
    @assert node_index >= 1 && node_index <= n_nodes "node_index out of bounds"
    @assert length(node_names) == n_nodes "node_names length must match matrix size"
    
    # Extract regulators for this node (column node_index)
    weights = interaction_matrix[:, node_index]
    
    # Identify positive and negative regulators (excluding self)
    pos_indices = Int[]
    neg_indices = Int[]
    pos_weights = Float64[]
    neg_weights = Float64[]
    
    for i in 1:n_nodes
        if i == node_index
            continue  # Skip self for now
        end
        if weights[i] > 0
            push!(pos_indices, i)
            push!(pos_weights, weights[i])
        elseif weights[i] < 0
            push!(neg_indices, i)
            push!(neg_weights, weights[i])
        end
    end
    
    P = length(pos_indices)
    N = length(neg_indices)
    SR = weights[node_index]  # self-regulation
    
    # Determine theta
    if theta === nothing
        if theta_mode == :half_difference
            theta = (P - N) / 2.0
        elseif theta_mode == :zero
            theta = 0.0
        else
            error("Must provide theta or use valid theta_mode")
        end
    elseif isa(theta, Function)
        theta = theta(P, N, SR)
    end
    
    # Get weights (assume uniform for now, but could be extracted from matrix)
    w_pos = P > 0 ? pos_weights[1] : 1.0
    w_neg = N > 0 ? -1*abs(neg_weights[1]) : -1.0
    
    # Calculate threshold conditions
    conditions = Boolean.calculate_threshold_conditions_general(
        P, N, SR; theta=theta,
        w_pos=w_pos, w_neg=w_neg, state_values=state_values
    )
    
    # Generate DNF clauses
    dnf_clauses = Vector{String}[]
    target_node_name = node_names[node_index]
    
    # Filter conditions where output is s_high (we want output = 1 or active state)
    active_conditions = filter(c -> (c.si_condition == :none && c.output == s_high) ||
                                     (c.si_condition == :positive && c.si_high_output == s_high) ||
                                     (c.si_condition == :negative && c.si_low_output == s_high),
                              conditions)
    
    for cond in active_conditions
        P_plus = cond.P_plus
        N_plus = cond.N_plus
        
        # Generate all combinations of choosing P_plus from P activators
        pos_combinations = P > 0 ? collect(combinations(1:P, P_plus)) : [[]]
        
        # Generate all combinations of choosing N_minus from N inhibitors (to be inactive)
        neg_combinations = N > 0 ? collect(combinations(1:N, N_plus)) : [[]]
        
        # For each combination, create an AND clause
        for pos_combo in pos_combinations
            for neg_combo in neg_combinations
                clause = String[]
                
                # Add positive regulators
                for i in 1:P
                    if i in pos_combo
                        # This activator should be ON (active)
                        push!(clause, node_names[pos_indices[i]])
                    else
                        # This activator should be OFF (inactive)
                        push!(clause, "~" * node_names[pos_indices[i]])
                    end
                end
                
                # Add negative regulators
                for i in 1:N
                    if i in neg_combo
                        # This inhibitor should be OFF (inactive) - positive literal
                        push!(clause, node_names[neg_indices[i]])
                    else
                        # This inhibitor should be ON (active) - negative literal
                        push!(clause, "~" * node_names[neg_indices[i]])
                    end
                end
                
                # Add self-regulation condition if needed
                if cond.si_condition == :positive
                    # Need current state to be high
                    push!(clause, target_node_name)
                elseif cond.si_condition == :negative
                    # Need current state to be low
                    push!(clause, "~" * target_node_name)
                end
                
                # Only add non-empty clauses
                if !isempty(clause)
                    push!(dnf_clauses, clause)
                end
            end
        end
    end
    
    # Get regulator names
    regulator_names = vcat(
        node_names[pos_indices],
        node_names[neg_indices]
    )
    
    # Add self if it appears in any clause
    if SR != 0 && any(c -> c.si_condition != :none, active_conditions)
        if !(target_node_name in regulator_names)
            push!(regulator_names, target_node_name)
        end
    end
    
    return dnf_clauses, regulator_names
end


"""
    interaction_matrix_to_dnf_all_nodes(
        interaction_matrix::Matrix,
        node_names::Vector{String};
        kwargs...
    ) -> Dict{String, Tuple{Vector{Vector{String}}, Vector{String}}}

Convert interaction matrix to DNF for all nodes in the network.

Returns a dictionary mapping node_name => (dnf, regulators)
"""
function interaction_matrix_to_dnf_all_nodes(
    interaction_matrix::Matrix,
    node_names::Vector{String};
    kwargs...
)
    n_nodes = size(interaction_matrix, 1)
    result = Dict{String, Tuple{Vector{Vector{String}}, Vector{String}}}()
    
    for i in 1:n_nodes
        dnf, regs = interaction_matrix_to_dnf(
            interaction_matrix, node_names, i; kwargs...
        )
        result[node_names[i]] = (dnf, regs)
    end
    
    return result
end


"""
    dnf_to_boolean_string(dnf::Vector{Vector{String}}) -> String

Convert DNF to readable Boolean expression string.

# Example
```julia
dnf = [["A", "B"], ["~A", "C"]]
expr = dnf_to_boolean_string(dnf)  # Returns "(A AND B) OR (~A AND C)"
```
"""
function dnf_to_boolean_string(dnf::Vector{Vector{String}})
    if isempty(dnf)
        return "FALSE"
    end
    
    terms = String[]
    for clause in dnf
        if isempty(clause)
            # push!(terms, "TRUE")
        else
            term = join([replace(literal, "~" => " NOT ") for literal in clause], " AND ")
            push!(terms, "($term)")
        end
    end
    
    return join(terms, " OR ")
end

function topo_to_boolean_rules(topo_file::String; kwargs...)
    J, nodes = topo2interaction(topo_file)
    rule = String[]
    for i in eachindex(nodes)
        dnf, regulators = interaction_matrix_to_dnf(J, nodes, i; kwargs...)
        if isempty(dnf)
            continue
        end
        dnf = simplify_dnf(dnf)
        target = nodes[i]
        boolean_string = dnf_to_boolean_string(dnf)
        push!(rule, target * " = " * boolean_string)
    end
    file_name = replace(topo_file, ".topo" => "_ising.txt")
    io = open(file_name, "w")
    for i in rule
        println(io, i)
    end
    close(io);
end