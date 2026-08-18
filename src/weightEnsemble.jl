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
