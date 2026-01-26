# Run benchmarks for different parameters
using BenchmarkTools
include("../src/Boolean.jl")
using .Boolean
using DataFrames
using CSV
results = Dict{Tuple{String, Symbol}, BenchmarkTools.Trial}()
rules_files = ["EMT_Switch.txt", "phase_switch.txt", "Phase_Switch_NEW.txt"]
methods = [:trajectory, :dnf]
for rules_file in rules_files
    for mtd in methods
        results[(rules_file, mtd)] = @benchmark simulate_network_logical($rules_file; 
        folder = "/Users/kishorehari/Desktop/Boolean.jl/data", 
        n_initial_conditions=10000, exact = false, method = $mtd) samples=100 seconds=300 evals=1
    end
end

# Create comparison table
comparison = DataFrame(
    Network = String[],
    n_nodes = Int[],
    n_edges = Int[],
    Method = Symbol[],
    median_ms = Float64[],
    memory_MiB = Float64[],
    allocs = Int[]
)
cd("/Users/kishorehari/Desktop/Boolean.jl/data")
for ((Network, Method), r) in results
    # Network, Method = n
    nodes = readlines(replace(Network, ".txt" => "_nodes.txt"))
    topoFile = readlines(replace(Network, ".txt" => ".topo"))
    println(r)
    push!(comparison, (
        Network = Network,
        n_nodes = length(nodes),
        n_edges = length(topoFile),
        Method = Method,
        median_ms = median(r.times) / 1e6,
        memory_MiB = r.memory / 1024^2,
        allocs = r.allocs
    ))
end
cd("..")
println(comparison)
CSV.write("scaling_benchmark.csv", comparison)