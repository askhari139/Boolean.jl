#=
Author : Kishore Hari
function name : syncUpdate
Description : simulate a network using ising formalism for multiple random initial conditions
Inputs :
    update_matrix : 2-D array; adjecency matrix of the network. each ij element depicts the edge from node i to node j.
    nInit : integer; number of initial conditions to use for simulation
    nIter : integer; maximum number of time steps for simulation
    
Active outputs :
    state_df : DataFrame; 3 columns - [init, fin, flag]:
        init : string; Initial state
        fin : string; Final state
        flag : integer; 1 - fin is a steady state, 0 - not
    Nodes : string array; List of nodes with the order in which state is listed
Passive outputs :
    The function writes states_df to a csv file if csv = true in input
=#
function asyncUpdate(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, stateRep::Int)
    n_nodes = size(update_matrix,1)
    stateVec = ifelse(stateRep == 0, [0,1], [-1,1])
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_matrix2 = sparse(update_matrix2')
    updFunc = ifelse(stateRep == 0, zeroConv, signVec)
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(zeroConv(state), "_")
        flag = 0
        time = 0
        j = 1
        while j <= nIter
            time = time + 1
            s1 = updFunc(update_matrix2*state)
            if iszero(j%10) # check after every ten steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            while s1[u] == state[u]
                
                j = j + 1
                time = time + 1
                if j > nIter
                    break
                end
                u = uList[j]
            end
            if flag == 1
                break
            end
            state[u] = s1[u]
        end
        if stateRep == 0
            fr = frustration(state, findnz(sparse(update_matrix)); negConv = true)
        else
            fr = frustration(state, findnz(sparse(update_matrix)))
        end
        fin = join(zeroConv(state), "_")
        push!(frustVec, fr)
        push!(initVec, init)
        push!(finVec, fin)
        push!(flagVec, flag)
        push!(timeVec, time)       
        # push!(states_df, (init, fin, flag))
    end
    states_df = DataFrame(init=initVec, 
            fin = finVec, flag = flagVec)
    frust_df = DataFrame(fin = finVec, 
        frust = frustVec)
    frust_df = unique(frust_df, :fin)
    timeData = DataFrame(fin = finVec, time = timeVec)
    timeData = groupby(timeData, :fin)
    timeData = combine(timeData, :time => avg, renamecols = false)
    frust_df = innerjoin(frust_df, timeData, on = :fin)
    # print(frust_df)
    return states_df, frust_df
end