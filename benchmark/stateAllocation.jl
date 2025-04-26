using BenchmarkTools
function generate_dense_matrix(n_nodes::Int, density::Float64, values = [-1, 1])
    A = zeros(Int8, n_nodes, n_nodes)
    for i in 1:n_nodes, j in 1:n_nodes
        if rand() < density
            A[i, j] = rand(values)
        end
    end
    return A
end
# Parameters
n_nodes = 50
nInit = 1000

# Dense update matrix
A = generate_dense_matrix(n_nodes, 0.2)

# Simulate stateList: n_nodes × nInit
stateVec = [-1, 1]
stateList = [rand(stateVec, n_nodes) for _ in 1:nInit]

# Regular signVec
signVec(x) = sign.(x)

# In-place signVec
function signVec!(out::Vector{Int}, x::Vector{Int})
    @inbounds for i in eachindex(x)
        out[i] = sign(x[i])
    end
    return out
end

# Preallocated buffers
state = Vector{Int}(undef, n_nodes)
s1 = similar(state)

println("=== Benchmark 1: Copy vs Allocate for `state` ===")
@btime state2 = $stateList[1];  # allocates
@btime copyto!($state, $stateList[1]);  # no alloc

println("\n=== Benchmark 2: Dense matvec + signVec (allocating) ===")
@btime begin
    s = $A * $state
    out = signVec(s)
end

println("\n=== Benchmark 3: In-place matvec + signVec! (non-allocating) ===")
@btime begin
    s1 = $A * $state     # in-place matrix-vector multiplication
    signVec!($s1, $s1)        # overwrite s1 in-place
end