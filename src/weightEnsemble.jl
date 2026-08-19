"""
    buildWeightStationaryDist(randVec0; noiseType="Multiplicative", noise=0.01,
                               nSteps=50_000, burnIn=10_000)

Simulate the per-edge weight random walk used as the RACIPE-"lambda" analog for
a weighted-network variant. `randVec0` is a base variant's per-edge magnitude
vector in [0,1] (the same convention as `edgeWeightPert`'s `_edgeWeightPert.dat`
columns / `bmodel_reps`'s `randVec` argument, where the *sign* of each edge is
carried separately by the underlying `update_matrix` and only the magnitude is
perturbed here).

Mirrors `RACIPEdata/src/simulate_racipe_resampled.jl`'s `build_lambda_dist`:
same clamped noise formulas as `asyncRandCont` (`async_update.jl`), just run as
a standalone per-edge walk (no Boolean state dynamics) since only the
stationary distribution of the weight itself is needed. Because magnitudes
live in [0,1] with no sign ambiguity, "Multiplicative" always targets 1
(the analog of `target = update_matrix_original[l] > 0 ? 1 : -1` when the sign
has already been factored out).

Returns a `(nSteps - burnIn) x nEdges` matrix of post-burn-in stationary
samples, one column per edge (matching `randVec0`'s edge order).
"""
function buildWeightStationaryDist(randVec0::Vector{Float64};
    noiseType::String="Multiplicative", noise::Float64=0.01,
    nSteps::Int=50_000, burnIn::Int=10_000)

    nEdges = length(randVec0)
    w = copy(randVec0)
    w0 = copy(randVec0)
    samples = Matrix{Float64}(undef, nSteps - burnIn, nEdges)

    for step in 1:nSteps
        xi = randn(nEdges) .* noise
        if noiseType == "Additive"
            w = w .+ xi
        elseif noiseType == "Multiplicative"
            w = w .+ (1.0 .- w) .* xi
        elseif noiseType == "Fluctuating"
            w = w0 .+ xi
        else
            error("Unknown noiseType '$noiseType'")
        end
        w = clamp.(w, 0.0, 1.0)
        if step > burnIn
            samples[step - burnIn, :] = w
        end
    end
    return samples
end

"""
    resampleRandVecs(stationarySamples; nResample=100)

Draw `nResample` new weight vectors, each edge drawn *independently* from its
own column of `stationarySamples` (mirrors `RACIPEdata`'s
`sample_param_vector`, which replaces each lambda with an independent draw
from that lambda's own stationary distribution — no cross-edge correlation is
preserved).
"""
function resampleRandVecs(stationarySamples::Matrix{Float64}; nResample::Int=100)
    nSamples, nEdges = size(stationarySamples)
    resampled = Matrix{Float64}(undef, nResample, nEdges)
    for k in 1:nEdges
        idx = rand(1:nSamples, nResample)
        resampled[:, k] = stationarySamples[idx, k]
    end
    return resampled
end

"""
    detResamplePert(topoFile, baseRandVecs; noiseType="Multiplicative",
                     noiseLevels=[0.001,0.005,0.01,0.02,0.05], nResample=100,
                     nInit=10_000, nIter=1_000, reps=3, rootPrefix="")

Deterministic ("det") arm counterpart to running `contWeightPert` (the
"stoch" arm) on the same base weighted-network variants. For each row of
`baseRandVecs` (one base variant's per-edge magnitude vector) and each noise
level: build that variant's per-edge stationary weight distribution
(`buildWeightStationaryDist`), draw `nResample` new variants from it
(`resampleRandVecs`), and run `bmodel_reps` (deterministic, no dynamics noise)
on each — writing one `*_finFlagFreq.csv` per resample via the existing
`randSim=true, randVec=...` path (same mechanism `edgeWeightPert` already
uses).
"""
function detResamplePert(topoFile::String, baseRandVecs::Matrix{Float64};
    noiseType::String="Multiplicative",
    noiseLevels::Vector{Float64}=[0.001, 0.005, 0.01, 0.02, 0.05],
    nResample::Int=100, nInit::Int=10_000, nIter::Int=1_000, reps::Int=3,
    rootPrefix::String="")

    nBase = size(baseRandVecs, 1)
    for b in 1:nBase
        for noise in noiseLevels
            dist = buildWeightStationaryDist(baseRandVecs[b, :]; noiseType=noiseType, noise=noise)
            resampled = resampleRandVecs(dist; nResample=nResample)
            p = Progress(nResample)
            Threads.@threads for r in 1:nResample
                root = "$(rootPrefix)base$(b)_$(noiseType)_$(noise)_$(r)"
                bmodel_reps(topoFile; nInit=nInit, nIter=nIter, randSim=true,
                    randVec=resampled[r, :], root=root, reps=reps)
                next!(p)
            end
            finish!(p)
        end
    end
end

"""
    writeWeightedTopo(topoFile, randVec, outFile)

Write a copy of `topoFile` with `Type` replaced by signed weights
(`sign(Type) * randVec[i]`, using the same nonzero-edge ordering as
`edgeWeightPert`/`bmodel_reps`'s `randVec` convention: `enumerate(findall(update_matrix .!= 0))`
over the matrix built by `topo2interaction`). Lets a base weighted-network
variant (otherwise only represented as a `randVec`, for the deterministic
`bmodel_reps` path) be fed to `contWeightPert` for the stochastic
(dynamics-noise) arm on the *same* base variant used in the deterministic arm.
"""
function writeWeightedTopo(topoFile::String, randVec::Vector{Float64}, outFile::String)
    df = DataFrame(CSV.File(topoFile))
    dropmissing!(df)
    update_matrix, Nodes = topo2interaction(topoFile)
    nzId = collect(enumerate(findall(update_matrix .!= 0)))
    idxMap = Dict((Nodes[j[1]], Nodes[j[2]]) => i for (i, j) in nzId)
    newType = Float64[]
    for row in eachrow(df)
        i = idxMap[(row[1], row[2])]
        push!(newType, sign(row[3]) * randVec[i])
    end
    outDf = copy(df)
    outDf[!, 3] = newType
    CSV.write(outFile, outDf; delim=" ")
end

@inline function teamScore(state::Vector{Float64}, idx::Vector{Int}, digits::Int)
    s = 0.0
    @inbounds for i in idx
        s += state[i] > 0 ? 1.0 : 0.0
    end
    return round(s / length(idx), digits=digits)
end

"""
    runTrajectoryScore!(stateMatrix, i, tid, localDict, localNextIDs, initState,
                        update_matrix, update_matrix_original, nzId, n_nodes,
                        nIter, nRecorded, frequency, noise, noiseType, stateRep,
                        team1Idx, team2Idx, scoreDigits, saveAt)

The per-trajectory body of `asyncRandContScore`'s `Threads.@threads for i in
1:nInit` loop, pulled out into its own function with fully concrete argument
types. Running this logic directly inside the `Threads.@threads` closure (as
`asyncRandCont` does) measured 800M+ allocations / ~20GiB for a mere
100x5000-step run despite the hot loop itself doing simple scalar/small-array
work -- a classic Julia performance pitfall where a closure capturing many
outer-scope variables (some behind a runtime `Union` such as `initialStates`,
which is `Vector{String}` or `Vector{Vector{Float64}}` depending on
`steadyStates`) can't be fully type-inferred, so ordinary operations get
boxed. Passing everything in as concretely-typed arguments to a dedicated
function lets Julia specialize and compile it properly.
"""
function runTrajectoryScore!(stateMatrix::Matrix{Int}, i::Int, tid::Int,
    localDict::Dict{Tuple{Float64,Float64},Int}, localNextIDs::Vector{Int},
    initState::Vector{Float64},
    update_matrix::Matrix{Float64}, update_matrix_original::Matrix{Float64},
    nzId::Vector{Tuple{Int,CartesianIndex{2}}},
    n_nodes::Int, nIter::Int, nRecorded::Int,
    frequency::Int, noise::Float64, noiseType::String, stateRep::Int,
    team1Idx::Vector{Int}, team2Idx::Vector{Int}, scoreDigits::Int, saveAt::Int)

    state = copy(initState)
    localUpdateMatrix = copy(update_matrix)
    uList = rand(1:n_nodes, nIter)

    scoreHistory = Matrix{Float64}(undef, nRecorded, 2)
    scoreHistory[1, 1] = teamScore(state, team1Idx, scoreDigits)
    scoreHistory[1, 2] = teamScore(state, team2Idx, scoreDigits)
    rIdx = 1

    for j in 2:nIter
        if iszero(j % frequency)
            randVec = randn(length(nzId)) * noise
            for (k, l) in nzId
                if noiseType == "Additive"
                    rVal = localUpdateMatrix[l] + randVec[k]
                elseif noiseType == "Multiplicative"
                    target = update_matrix_original[l] > 0 ? 1 : -1
                    rVal = localUpdateMatrix[l] + (target - localUpdateMatrix[l]) * randVec[k]
                elseif noiseType == "Fluctuating"
                    rVal = update_matrix_original[l] + randVec[k]
                else
                    rVal = localUpdateMatrix[l] + randVec[k]
                end
                rVal = localUpdateMatrix[l] > 0 ? min(max(rVal, 0), 1) : min(max(rVal, -1), 0)
                localUpdateMatrix[l] = rVal
            end
        end

        u = uList[j]
        raw = 0.0
        @inbounds for k in 1:n_nodes
            raw += localUpdateMatrix[u, k] * state[k]
        end
        state[u] = scalarUpdate(raw, stateRep, state[u])
        if iszero((j - 1) % saveAt)
            rIdx += 1
            scoreHistory[rIdx, 1] = teamScore(state, team1Idx, scoreDigits)
            scoreHistory[rIdx, 2] = teamScore(state, team2Idx, scoreDigits)
        end
    end

    for j in 1:nRecorded
        key = (scoreHistory[j, 1], scoreHistory[j, 2])
        if !haskey(localDict, key)
            localDict[key] = localNextIDs[tid]
            localNextIDs[tid] += 1
        end
        stateMatrix[i, j] = localDict[key]
    end
    return nothing
end

"""
    asyncRandContScore(update_matrix, nInit, nIter, stateRep, team1Idx, team2Idx;
                        randVec=[0.0], noise=0.01, frequency=1, steadyStates=true,
                        topN=10, noiseType="Additive", noiseMode="All", scoreDigits=1,
                        saveAt=1)

Same noisy asynchronous Boolean dynamics as `asyncRandCont`, but tracks MRT at
the level of the `(T1, T2)` team-score bin instead of the raw per-node Boolean
state. `T1`/`T2` are each `round(mean(zeroConv(state[teamIdx])), digits=scoreDigits)`
(`zeroConv` maps the internal -1/1 or 0/1 state to 0/1, matching this project's
existing `eScore`/`mScore` convention).

Raw-state tracking is the wrong granularity for this pipeline: with random
initial conditions and no attractor-seeding, a single 1000x10000-trajectory
run visited 558,867 distinct raw states (mostly transients visited once),
each paying for a dictionary insert + eventual string conversion, even though
MRT is only ever consumed at the (T1,T2) bin level (at most
`(1/10^-scoreDigits * 2 + 1)^2` bins, a few hundred at most here). Scoring
inline instead of after the fact collapses the dictionary accordingly, and
also avoids storing the full per-step state history (`scoreHistory` is
`nRecorded x 2`, not `nIter x n_nodes`).

This is a separate function rather than an option on `asyncRandCont` so that
existing raw-state callers are completely unaffected (no dict-type or
branching overhead added to that hot loop).

`frequency` is *not* a sampling-rate knob -- it's how often (in async steps)
the edge weights themselves get re-perturbed (the existing `asyncRandCont`
noise mechanism, here doubling as the RACIPE-"lambda"-drift-rate analog, i.e.
what a user thinking in ODE terms would call the network's `DT`). At
`frequency=1` the network is re-perturbed on every single-node async update,
leaving the rest of the network no time to relax under a fixed set of
weights before they shift again; `frequency=10` or `100` lets `~frequency`
async updates happen (roughly `frequency / n_nodes` full sweeps of the
network) between re-perturbations.

`saveAt` is the actual sampling-rate knob (like an ODE solver's `saveat`):
only every `saveAt`-th step (plus the initial condition) is scored and
tallied into the MRT dictionary, rather than all `nIter` steps -- both
because consecutive raw steps are highly autocorrelated (only one node
differs) and because it cuts the number of `teamScore` calls and dictionary
inserts by the same factor. `stateMatrix` therefore has
`1 + fld(nIter - 1, saveAt)` columns, not `nIter`.

Returns `(stateMatrix, sListUnique)` exactly like `asyncRandCont`, except
`sListUnique` entries are `"T1_T2"` labels (e.g. `"0.3_-0.5"`) instead of
per-node bit-strings, and `stateMatrix` has one column per *recorded* step
rather than per raw step.
"""
function asyncRandContScore(update_matrix::Union{Array{Int,2}, Array{Float64,2}},
    nInit::Int, nIter::Int, stateRep::Int, team1Idx::Vector{Int}, team2Idx::Vector{Int};
    randVec::Array{Float64,1} = [0.0], noise::Float64 = 0.01,
    frequency::Int = 1, steadyStates::Bool = true, topN::Int = 10,
    noiseType::String = "Additive", noiseMode::String = "All", scoreDigits::Int = 1,
    saveAt::Int = 1)

    n_nodes = size(update_matrix, 1)

    if steadyStates
        updm = Int.(sign.(update_matrix))
        states_df, frust_df = asyncUpdate(updm, 10000, 1000, stateRep, false, Int[], Int[], Int[])
        if topN != 0
            freq_table = combine(groupby(states_df, :fin), nrow => :Count)
            freq_table = sort(freq_table, :Count, rev=true)
            states = freq_table[1:min(topN, nrow(freq_table)), :fin]
        else
            states = states_df[:, :fin]
        end
        initialStates = rand(states, nInit)
    else
        stateVec = ifelse(stateRep == 0, [0.0, 1.0], [-1.0, 1.0])
        initialStates = [rand(stateVec, n_nodes) for _ in 1:nInit]
    end

    update_matrix = update_matrix'
    nzIdOrig = enumerate(findall(update_matrix .!= 0))
    nzId = typeof(collect(nzIdOrig))()
    noiseMode = noiseMode in ["All", "Activation", "Inhibition"] ? noiseMode : "All"
    if noiseMode == "All"
        nzId = collect(nzIdOrig)
    else
        for (k, l) in nzIdOrig
            if noiseMode == "Activation"
                update_matrix[l] > 0 && push!(nzId, (k, l))
            else
                update_matrix[l] < 0 && push!(nzId, (k, l))
            end
        end
    end

    update_matrix = Matrix{Float64}(float(update_matrix))
    update_matrix_original = copy(update_matrix)

    nRecorded = 1 + fld(nIter - 1, saveAt)
    stateMatrix = zeros(Int, nInit, nRecorded)

    nthreads = Threads.nthreads()
    localDicts = [Dict{Tuple{Float64,Float64}, Int}() for _ in 1:nthreads]
    localNextIDs = ones(Int, nthreads)
    trajectoryThreads = zeros(Int, nInit)

    Threads.@threads for i in 1:nInit
        tid = Threads.threadid()
        trajectoryThreads[i] = tid
        initState::Vector{Float64} = steadyStates ? parse.(Float64, split(initialStates[i], "_")) : initialStates[i]
        runTrajectoryScore!(stateMatrix, i, tid, localDicts[tid], localNextIDs,
            initState, update_matrix, update_matrix_original, nzId, n_nodes, nIter, nRecorded,
            frequency, noise, noiseType, stateRep, team1Idx, team2Idx, scoreDigits, saveAt)
    end

    globalDict = Dict{Tuple{Float64,Float64}, Int}()
    globalNextID = 1
    for tid in 1:nthreads
        for (key, local_id) in localDicts[tid]
            if !haskey(globalDict, key)
                globalDict[key] = globalNextID
                globalNextID += 1
            end
        end
    end

    Threads.@threads for i in 1:nInit
        original_tid = trajectoryThreads[i]
        localDict = localDicts[original_tid]
        reverseMap = Dict(local_id => key for (key, local_id) in localDict)
        for j in 1:nRecorded
            local_id = stateMatrix[i, j]
            key = reverseMap[local_id]
            stateMatrix[i, j] = globalDict[key]
        end
    end

    sListUnique = Vector{String}(undef, length(globalDict))
    for (key, id) in globalDict
        sListUnique[id] = string(key[1], "_", key[2])
    end

    return stateMatrix, sListUnique
end

"""
    contWeightPertScore(topoFile, team1, team2; nInit=1000, nIter=100000,
                         mode="Async", stateRep=-1, noise=0.01,
                         steadyStates=true, topN=10, frequency=1,
                         noiseType="Additive", noiseMode="All", scoreDigits=1,
                         saveAt=1)

Same outputs as `contWeightPert` (`*_contWeightPert.csv`,
`*_contWeightPert_states.csv`, `*_contWeightPert_trajectories.csv`,
`*_contWeightPert_byInitialState.csv`, with a `_score` tag in the filename to
keep them distinguishable from `contWeightPert`'s raw-state outputs), but MRT
is tracked at the `(T1, T2)` team-score bin level via `asyncRandContScore`
instead of the raw Boolean state -- the `State`/`InitialState` columns
contain `"T1_T2"` labels instead of per-node bit-strings. `team1`/`team2` are
node *names* (as in a `.teams` file), resolved to indices via
`topo2interaction`'s node ordering.

`frequency` (the network re-perturbation rate, i.e. "DT") and `saveAt` (the
MRT-sampling rate) both change how many columns `asyncRandContScore` actually
returns (`1 + fld(nIter - 1, saveAt)`, not `nIter`), and both are folded into
the output filename (`_DT<frequency>_SA<saveAt>_score`) so sweeps over either
don't collide.
"""
function contWeightPertScore(topoFile::String, team1::Vector{String}, team2::Vector{String};
    nInit::Int64=1000, nIter::Int64=100000, mode::String="Async", stateRep::Int64=-1,
    noise::Float64=0.01, steadyStates::Bool=true, topN::Int=10,
    frequency::Int=1, noiseType::String="Additive", noiseMode::String="All",
    scoreDigits::Int=1, saveAt::Int=1)

    updMat, nodes = topo2interaction(topoFile)
    team1Idx = [findfirst(==(n), nodes) for n in team1]
    team2Idx = [findfirst(==(n), nodes) for n in team2]

    y1 = @elapsed stateMat, states = asyncRandContScore(updMat, nInit, nIter, stateRep,
        team1Idx, team2Idx; noise=noise, steadyStates=steadyStates, topN=topN,
        frequency=frequency, noiseType=noiseType, noiseMode=noiseMode, scoreDigits=scoreDigits,
        saveAt=saveAt)
    println("Async rand (score) took $y1")
    nStates = length(states)
    nCols = size(stateMat, 2)

    stateFreqs = zeros(Float64, nInit, nStates)
    Threads.@threads for i in 1:nInit
        counts = Dict{Int, Int}()
        for j in 1:nCols
            stateID = stateMat[i, j]
            counts[stateID] = get(counts, stateID, 0) + 1
        end
        for (stateID, count) in counts
            stateFreqs[i, stateID] = count / nCols
        end
    end

    MRT = vec(mean(stateFreqs, dims=1))
    MRTsd = vec(std(stateFreqs, dims=1))

    switchingEvents = zeros(Int, nInit)
    for i in 1:nInit
        switchingEvents[i] = sum(stateMat[i, 2:end] .!= stateMat[i, 1:end-1])
    end

    baseName = replace(topoFile, ".topo" => "")
    baseName = baseName*"_"*noiseType*"_"*(if noiseMode != "All" noiseMode*"_" else "" end)*string(noise)*
        "_DT"*string(frequency)*"_SA"*string(saveAt)*"_score"

    trajectoryDF = DataFrame(stateMat, :auto)
    rename!(trajectoryDF, ["t$i" for i in 1:nCols])
    CSV.write("$(baseName)_contWeightPert.csv", trajectoryDF)

    stateStatsDF = DataFrame(ID = 1:nStates, State = states, MRT = MRT, MRTsd = MRTsd)
    sort!(stateStatsDF, :MRT, rev=true)
    CSV.write("$(baseName)_contWeightPert_states.csv", stateStatsDF)

    initialStateIDs = stateMat[:, 1]
    initialStates = states[initialStateIDs]

    trajectoryStatsDF = DataFrame(
        TrajectoryID = 1:nInit,
        InitialStateID = initialStateIDs,
        InitialState = initialStates,
        SwitchingEvents = switchingEvents
    )
    CSV.write("$(baseName)_contWeightPert_trajectories.csv", trajectoryStatsDF)

    if steadyStates
        initialConditionStats = combine(groupby(trajectoryStatsDF, [:InitialStateID, :InitialState])) do df
            trajIDs = df.TrajectoryID
            stateMRTs = [mean(stateFreqs[trajIDs, stateID]) for stateID in 1:nStates]
            topIndices = sortperm(stateMRTs, rev=true)[1:min(5, nStates)]
            DataFrame(
                nTrajectories = nrow(df),
                MeanSwitchingEvents = mean(df.SwitchingEvents),
                StdSwitchingEvents = std(df.SwitchingEvents),
                MedianSwitchingEvents = median(df.SwitchingEvents),
                TopStates = join(states[topIndices], "; "),
                TopStatesMRT = join(round.(stateMRTs[topIndices], digits=4), "; ")
            )
        end
        CSV.write("$(baseName)_contWeightPert_byInitialState.csv", initialConditionStats)
    end

    return stateStatsDF, trajectoryStatsDF
end
