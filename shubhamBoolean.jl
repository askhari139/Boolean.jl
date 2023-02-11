function shubhamBoolean(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int)
    n_nodes = size(update_matrix,1)
    stateVec = Int[-1,1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_matrix2 = sparse(update_matrix2')
    for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        flag = 0
        for j in 1:nIter
            s1 = sign.(update_matrix2*state)
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = frustration(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        push!(frustVec, fr)
        push!(initVec, init)
        push!(finVec, fin)
        push!(flagVec, flag)       
        # push!(states_df, (init, fin, flag))
    end
    states_df = DataFrame(init=initVec, 
            fin = finVec, flag = flagVec)
    frust_df = DataFrame(fin = finVec, 
        frust = frustVec)
    return states_df, unique(frust_df, :fin)
end
