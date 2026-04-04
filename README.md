# Boolean.jl

A Julia package for simulating Boolean and multi-level gene regulatory networks. Boolean.jl supports Ising-type, non-Ising (threshold-based), logical-rule, and continuous-weight formalisms, and provides tools for attractor identification, perturbation analysis, stochastic trajectory simulation, and State Transition Graph (STG) construction.

---

## Table of Contents

1. [Installation](#installation)
2. [Input File Formats](#input-file-formats)
3. [Core Functionalities](#core-functionalities)
4. [Quick Start](#quick-start)
5. [Detailed Usage](#detailed-usage)
   - [Ising / threshold-based simulation](#1-ising--threshold-based-simulation)
   - [Logical-rule simulation](#2-logical-rule-simulation)
   - [Stochastic (continuous-weight) simulation](#3-stochastic-continuous-weight-simulation)
   - [Perturbation & coherence analysis](#4-perturbation--coherence-analysis)
   - [State Transition Graph (STG)](#5-state-transition-graph-stg)
   - [Network utilities](#6-network-utilities)
6. [Multi-threading](#multi-threading)
7. [Dependencies](#dependencies)
8. [Citation](#citation)

---

## Installation

### Option A — Add directly from GitHub

```julia
using Pkg
Pkg.add(url = "https://github.com/askhari139/Boolean.jl")
using Boolean
```


### Option B — Clone and load locally (recommended while the package is in active development)

```julia
# 1. Clone the repository
# git clone https://github.com/askhari139/Boolean.jl

# 2. From the Julia REPL, activate the environment and instantiate
using Pkg
Pkg.activate("/path/to/Boolean.jl")   # path to the cloned folder
Pkg.instantiate()                      # installs all dependencies

# 3. Load the package
using Boolean
```

Alternatively, use the provided loader script:

```julia
include("/path/to/Boolean.jl/load_local.jl")
using Boolean
```

### Compatibility

Tested on Julia **1.6, 1.8, 1.9, and 1.10**.

---

## Input File Formats

### Topology file (`.topo`)

A space-delimited edge list with a header row. Each subsequent row describes one directed edge.

```
Source Target Type
NodeA  NodeB  1
NodeC  NodeB  2
NodeB  NodeB  1
```

- **Type 1** — activation  
- **Type 2** — inhibition

### Boolean rules file (`.txt`)

One rule per line in the form `Target = Boolean expression`. Standard operators are `AND`, `OR`, `NOT`.

```
NodeA = NodeC OR NOT NodeD
NodeB = NodeA AND NOT NodeC
NodeC = NodeB
```

---

## Core Functionalities

| Functionality | Key function(s) |
|---|---|
| Ising / threshold-based attractor search | `bmodel_reps` |
| Logical-rule attractor search (sync / async) | `simulate_network_logical` |
| Stochastic continuous-weight simulation | `contWeightPert` |
| Discrete edge-weight perturbation scan | `edgeWeightPert` |
| Perturbation robustness / coherence | `coherence`, `coherenceAllNode` |
| Single-node knockout scan | `scanNodeTurnOff` |
| Stochastic trajectory generation | `simulate_async`, `simulate_multiple_states_to_df` |
| Logical trajectory generation | `simulate_multiple_states_to_df_logical` |
| State Transition Graph construction | `topo_to_stg`, `bool_to_stg` |
| Topology ↔ Boolean rule conversion | `boolean_to_topo`, `topo_to_boolean_rules` |
| Extract node list from topo file | `getNodes` |

---

## Quick Start

```julia
using Boolean

# --- Ising simulation ---
bmodel_reps("myNetwork.topo"; nInit=10_000, nIter=1_000, mode="Async", stateRep=-1)
# Writes myNetwork_finFlagFreq.csv with attractor states and basin sizes

# --- Logical-rule simulation ---
simulate_network_logical("myNetwork.txt"; n_initial_conditions=10_000)
# Writes myNetwork_attractors.csv

# --- Stochastic continuous-weight simulation ---
contWeightPert("myNetwork.topo"; nInit=1_000, nIter=100_000, noise=0.01)
# Writes trajectory and state statistics CSVs
```

---

## Detailed Usage

### 1. Ising / threshold-based simulation

`bmodel_reps` is the main entry point for Ising-type simulations. It runs multiple replicates and aggregates attractor statistics.

```julia
bmodel_reps(
    "myNetwork.topo";
    nInit  = 10_000,   # initial conditions per replicate
    nIter  = 1_000,    # max update steps
    mode   = "Async",  # "Async" only (Sync is handled separately)
    stateRep = -1,     # -1 → {-1,1} Ising states; 0 → {0,1} non-Ising states
    reps   = 3,        # number of replicates
    types  = [0],      # 0 = equal weights; 1 = activation-dominant; 2 = inhibition-dominant
    write  = true,     # write output CSV
    getData = false    # return DataFrame instead of writing
)
```

**With knockdown / overexpression:**

```julia
# Node indices refer to the alphabetically sorted node list
bmodel_reps("myNetwork.topo"; kdNodes=[2], oeNodes=[5], nInit=10_000)
```

**Multi-level (Shubham) dynamics:**

```julia
bmodel_reps("myNetwork.topo"; shubham=true, nLevels=3, nInit=10_000)
```

**Output columns in `_finFlagFreq.csv`:**

| Column | Description |
|---|---|
| `states` | Attractor state string (`"0_1_0_1_..."`) |
| `flag` | 1 = converged, 0 = did not converge within `nIter` |
| `Avg0` | Mean basin size across replicates (for rule type 0) |
| `SD0` | Standard deviation of basin size |
| `frust0` | Network frustration score |
| `time` | Mean steps to convergence |

---

### 2. Logical-rule simulation

`simulate_network_logical` parses Boolean rules, builds truth tables (or DNF evaluators), and finds attractors under synchronous or asynchronous updates.

```julia
# Synchronous (default)
simulate_network_logical(
    "myNetwork.txt";
    update_mode           = "synchronous",
    n_initial_conditions  = 100_000,
    n_replicates          = 3,
    max_steps             = 1_000,
    exact                 = false   # true exhaustively covers all 2^N states (small networks only)
)
```

```julia
# Asynchronous
simulate_network_logical(
    "myNetwork.txt";
    update_mode           = "asynchronous",
    n_initial_conditions  = 100_000,
    max_steps             = 5_000
)
```

**Output columns in `_attractors.csv`:**

| Column | Description |
|---|---|
| `states` | Attractor cycle, states separated by `;` |
| `basin_size_mean` | Mean basin proportion across replicates |
| `basin_size_sd` | Standard deviation |
| `time_mean` | Mean convergence time |
| `frustration` | Network frustration score |
| `flag` | 1 = converged |

**Convert topology → Boolean rules:**

```julia
topo_to_boolean_rules("myNetwork.topo")
# Writes myNetwork_ising.txt with inferred Boolean rules
```

**Convert Boolean rules → topology:**

```julia
boolean_to_topo("myNetwork.txt")
# Writes myNetwork.topo
```

---

### 3. Stochastic (continuous-weight) simulation

`contWeightPert` adds Gaussian noise to edge weights at each step, enabling continuous exploration of the weight landscape.

```julia
contWeightPert(
    "myNetwork.topo";
    nInit        = 1_000,    # trajectories
    nIter        = 100_000,  # steps per trajectory
    noise        = 0.01,     # std dev of Gaussian noise on edge weights
    steadyStates = true,     # seed trajectories from known attractors
    topN         = 10,       # use the top N most-visited states as seeds
    stateRep     = -1
)
```

Outputs three CSVs:
- `*_contWeightPert.csv` — full trajectory matrix (trajectory × time)
- `*_contWeightPert_states.csv` — Mean Residence Time (MRT) per state
- `*_contWeightPert_trajectories.csv` — per-trajectory switching statistics

**Discrete edge-weight scan:**

```julia
edgeWeightPert(
    "myNetwork.topo";
    nPerts     = 10_000,   # number of random weight samples
    minWeight  = 0.0,
    maxWeight  = 1.0,
    nInit      = 10_000,
    nIter      = 1_000
)
```

---

### 4. Perturbation & coherence analysis

`coherence` measures how robustly each attractor recovers after single- or multi-node perturbations.

```julia
# Requires a prior bmodel_reps run (reads *_finFlagFreq.csv)
coherence(
    "myNetwork.topo";
    nPert   = 1,    # number of nodes to perturb simultaneously
    nInit   = 100,  # perturbation trials per attractor
    nSim    = 10,   # simulations per perturbed state
    nLevels = 1,    # 1 = Boolean {-1,1}; 2 = two-level multi-valued
    nIter   = 1_000
)
```

**Scan all perturbation sizes (1 to N):**

```julia
coherenceAllNode(
    "myNetwork.topo";
    nInit   = 100,
    nSim    = 10,
    nLevels = 1,
    nIter   = 1_000
)
```

**Single-node knockout scan:**

```julia
scanNodeTurnOff(
    "myNetwork.topo";
    nInit   = 10_000,
    nIter   = 1_000,
    nLevels = 2
)
# Writes *_scanNode_turnOff_finFlagFreq.csv
```

---

### 5. State Transition Graph (STG)

```julia
# Build STG from a topo file (Ising dynamics)
stg_ising  = topo_to_stg("myNetwork.topo"; mode=:ising)

# Build STG from Boolean rules
stg_logic  = bool_to_stg("myNetwork.txt")

# Extract attractors and basin sizes from an STG
attractor_dict = get_attractors_from_stg(stg_ising)
# Returns Dict{Vector{Vector{Int}}, Float64}  attractor_states => basin_fraction

# Compare logical vs Ising vs non-Ising dynamics
attractor_dict, jsd = compare_logical_ising("myNetwork.txt")
# jsd is a Dict of pairwise Jensen–Shannon divergences between the three formalisms
```

---

### 6. Network utilities

```julia
# Write node list to file (alphabetically sorted)
getNodes("myNetwork.topo")
# Writes myNetwork_nodes.txt

# Stochastic Ising trajectory (returns DataFrame)
df = simulate_multiple_states_to_df(
    "myNetwork.topo", 2;  # topo file, n (number of levels)
    nSim  = 10,
    steps = 1_000
)

# Logical trajectory (synchronous, with perturbations)
df = simulate_multiple_states_to_df_logical(
    "myNetwork.txt";
    nSim             = 10,
    perturb_nodes    = 2,
    perturb_interval = 200,
    max_steps        = 2_000
)
```

---

## Multi-threading

Most batch functions (`bmodel_reps`, `edgeWeightPert`, `scanNodeTurnOff`, `contWeightPert`) are parallelised with `Threads.@threads`. Launch Julia with multiple threads for significant speedups:

```bash
julia --threads 8 script.jl
# or set the environment variable
export JULIA_NUM_THREADS=8
julia script.jl
```

Check the number of active threads at runtime:

```julia
using Base.Threads
println(Threads.nthreads())
```

---

## Dependencies

Boolean.jl depends on the following registered Julia packages (installed automatically via `Pkg.instantiate`):

- `CSV`, `DataFrames` — tabular I/O and manipulation
- `SparseArrays`, `LinearAlgebra` — matrix operations
- `Random`, `Combinatorics` — stochastic simulation and combinatorics
- `StatsBase` — statistics utilities
- `ProgressMeter` — progress display
- `Pipe` — pipeline syntax
- `JSON` — serialisation
- `IterTools` — iterator utilities

---

## Citation

If you use Boolean.jl in your research, please cite the original Ising Boolean formalism:

> Font-Clos, F., Zapperi, S., & La Porta, C. A. M. (2018).
> Topography of epithelial–mesenchymal plasticity.
> *Proceedings of the National Academy of Sciences*, 115(23), 5902–5907.

And acknowledge this software repository:

> Hari, K. *Boolean.jl* — Julia code for Ising and Boolean network simulations.  
> https://github.com/askhari139/Boolean.jl
