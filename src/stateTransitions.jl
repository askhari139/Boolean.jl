# ============================
# Main asynchronous simulation function
# ============================
"""
    simulate_async!(s, J, n; steps=1000, rng, stateList, record_every,
                    perturb_nodes, perturb_mode, perturb_interval)

Simulate the multilevel Boolean network asynchronously **in-place** on integer-scaled state `s` (`-n..n`, no zeros).

# Arguments
- `s::Vector{Int}` : initial state vector (length = number of nodes), modified in-place
- `J::AbstractMatrix{Int}` : adjacency/sign matrix where `J[i,j]` is the influence of node `i` on node `j`
- `n::Int` : denominator (levels are `k/n`, k=-n..n, no zeros)

# Keyword Arguments
- `steps::Int=1000` : total number of asynchronous update steps
- `rng` : random number generator (default `Random.GLOBAL_RNG`)
- `stateList::Vector{Vector{Int}}=Vector{Vector{Int}}()` : optional history recording
- `record_every::Int=1` : record every `record_every` steps
- `perturb_nodes::Int=0` : number of nodes to perturb at each perturbation event
- `perturb_mode::Symbol=:mirror` : `:mirror` flips `k → -k`, `:full` flips to extreme `-n*sign(k)`
- `perturb_interval::Int=typemax(Int)` : steps between perturbations (`Inf` = no perturbations)

# Returns
`(s, stateList)` : final state vector and recorded history
"""
function simulate_async_transition!(
    s1::Vector{Int},
    s2::Vector{Int},
    J::AbstractMatrix{Int},
    n::Int;
    steps::Int = 1000,
    rng = Random.GLOBAL_RNG,
    record_every::Int = 1,
    perturb_levels::Int = 0,
    perturb_mode::Symbol = :mirror,
    perturb_interval::Int = typemax(Int)
)
    N = length(s1)
    @assert size(J,1) == N && size(J,2) == N "J must be square N×N matching s"

    diff_nodes = findall(v1 .!= v2)
    

    # Precompute denominators: d[i] = sum_j |J[j,i]|
    stateList = Vector{Vector{Int}}()
    push!(stateList, copy(s))
    d = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        ssum = 0
        for j in 1:N
            ssum += abs(J[j,i])
        end
        d[i] = ssum
    end

    # Main asynchronous loop
    @inbounds for t in 1:steps
        # pick one random node
        i = rand(rng, 1:N)
        di = d[i]
        if di > 0
            # compute input sum
            M = 0
            @inbounds for j in 1:N
                M += J[j,i] * s[j]
            end

            if M != 0
                if M < -(n-1)*di
                    newk = -n
                elseif M > (n-1)*di
                    newk = n
                else
                    if M < 0
                        k = (-M + di - 1) ÷ di
                        newk = -k
                    else
                        k = (M + di - 1) ÷ di
                        newk = k
                    end
                end
                # enforce non-zero state
                if newk == 0
                    newk = sign(M)  # ±1
                end
                s[i] = newk
            end
        end

        # Apply perturbation if scheduled
        if (perturb_interval > 0 && (t % perturb_interval == 0) || t == 10) && perturb_nodes > 0
            idxs = rand(rng, 1:N, perturb_nodes)  # sample with replacement
            for i in idxs
                if perturb_mode == :mirror
                    s[i] = -s[i]
                elseif perturb_mode == :full
                    s[i] = -n * sign(s[i])
                elseif perturb_mode == :next
                    s[i] = s[i] + rand([-1,1])
                else
                    error("Unknown perturbation mode: $perturb_mode")
                end
            end
        end

        # Record state
        if (t % record_every == 0)
            push!(stateList, copy(s))
        end
    end
    times = collect(0:(length(stateList)-1)) .* record_every
    return [s ./ n for s in stateList], times
end

# ============================
# Non-mutating wrapper
# ============================
"""
    simulate_async(J, n; s0, steps, rng, stateList, record_every,
                   perturb_nodes, perturb_mode, perturb_interval)

Non-mutating wrapper that copies `s0` or generates a random initial state, then calls `simulate_async!`.
"""
function simulate_async(
    J::AbstractMatrix{Int},
    n::Int;
    s0::Union{Nothing,AbstractVector{Int}} = nothing,
    steps::Int = 1000,
    rng = Random.GLOBAL_RNG,
    record_every::Int = 1,
    perturb_nodes::Int = 0,
    perturb_mode::Symbol = :mirror,
    perturb_interval::Int = typemax(Int)
)
    N = size(J,1)
    @assert size(J,2) == N "J must be square N×N"
    s_init = s0 === nothing ? random_state(N, n; rng=rng) : copy(collect(s0))
    return simulate_async!(
        s_init, J, n;
        steps=steps, rng=rng, record_every=record_every,
        perturb_nodes=perturb_nodes, perturb_mode=perturb_mode, perturb_interval=perturb_interval
    )
end
