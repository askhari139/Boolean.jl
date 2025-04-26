#=
Things I tried to improve performance:
1. Sparse matrix - doesn't make much difference. 
2. statelist generation - vector form of states perform much better than matrix form in indexing
3. static arrays - not much difference
4. replacing the dot operation in numChar with a loop - HUGE difference
5. Assigning size of the array in advance - not much difference

=#

using BenchmarkTools
topoFile = "EMT_RACIPE.topo"
update_matrix, nodes = Boolean.topo2interaction(topoFile)
n_nodes = size(update_matrix, 1)
updOriginal = copy(update_matrix)
update_matrix = float(update_matrix)
  # need it to calculate frustration
# dividing each column by the corresponding node indegree to ensure that the product is betweenn -1 and 1





s = Int64[]
@benchmark begin
    asyncUpdate(updOriginal, 100, 1000,-1,false,s,s,s,s)
end
@benchmark begin
    shubhamBoolean(updOriginal, 100, 1000,1,false,s,s,s)
end
@code_warntype asyncUpdate(updOriginal, 100, 1000,-1,false,s,s,s,s)
@code_warntype shubhamBoolean(updOriginal, 100, 1000,1,false,s,s,s)

for i in 1:n_nodes
    n = sum(update_matrix[:,i] .!= 0)        
    if n == 0
        n = 1
    end
    update_matrix[:, i] = update_matrix[:, i]/n
end

updOriginal = updOriginal'
update_matrix = update_matrix'
levels = [-1, -0.5, 0.5, 1]
sVecList = [levels for _ in 1:n_nodes]
vaibhav = [false for _ in 1:n_nodes]
levels = [levels for _ in 1:n_nodes]
# Assume these are initialized appropriately
@benchmark begin
    # single update for asyncUpdate
    state = rand([-1, 1], n_nodes)
    prod = updOriginal * state
    s1 = signVec(prod)
end

@benchmark begin
    # single update for shubhamBoolean
    state = rand([-1, 0.5, 0.5, 1], n_nodes)

    prod = update_matrix * state
    prod = stateChar(prod, state, levels, sVecList, vaibhav)
end

@benchmark begin
    # frustration calc
    state = rand([-1, 0.5, 0.5, 1], n_nodes)
    fr = shubhamFrust(state, findnz(sparse(update_matrix)))
end

@benchmark begin
    # string conversion
    state = rand([-1, 0.5, 0.5, 1], n_nodes)
    stateConvert(state, 2)
end
