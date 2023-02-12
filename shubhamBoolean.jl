function stateChar(state::AbstractArray)
    for i in 1:length(state)
        x = state[i]
        if x == 0
            y = copy(x)
        elseif x < -0.5
            y = -1.0
        elseif x < 0
            y = -0.5
        elseif x > 0.5
            y = 1.0
        else
            y = 0.5
        end
        state[i] = y
    end
    return state
end

function shubhamFrust(state::Array{Int,1}, 
    nonZeros::Tuple{Array{Int64,1},Array{Int64,1},Array{Float64,1}})
    frustration = 0
    nEdges = length(nonZeros[1])
    for (x,y,v) in zip(nonZeros...)
        s = state[x]*state[y]*v
        if s<0
            if state[x] == state[y]
                frustration = frustration + 1
            else
                frustration = frustration + 0.5
            end
        end
    end
    frustration = frustration/nEdges
    return frustration
end

function stateConvert(state)
    state = 2*state
    state = join(["'", join(replace(x -> x < 0 ? x+2 : x+1, state)), "'"])
    return state
end

function shubhamBoolean(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int)
    n_nodes = size(update_matrix,1)
    stateVec = [-1, -0.5, 0.5, 1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix = float(update_matrix)
    for i in 1:n_nodes
        n = length([(i,j) for j=1:n_nodes if u[i,j] == 0])        
        n = n_nodes - n
        update_matrix[i, :] = update_matrix[i,:]/n
    end
    # update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    # nzId = enumerate(findall(update_matrix.!=0))
    # minVal = minimum([abs(update_matrix[j]) for (i,j) in nzId])
    # update_matrix2 = update_matrix + Matrix(I, n_nodes, n_nodes)*(minVal/2)
    update_matrix2 = sparse(update_matrix2')
    for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        flag = 0
        for j in 1:nIter
            s1 = stateChar(update_matrix2*state)
            u = rand(1:n_nodes, 1)
            if iszero(j%2) # check after every two steps,hopefully reduce the time
                if s1 == state
                    flag = 1
                    break
                end
            end
            state[u] = s1[u]
        end
        fr = shubhamFrust(state, findnz(sparse(update_matrix)))
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
