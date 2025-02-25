# Function that takes in an update matrix and returns a set of N functions corresponding to each node of the network

function updateFuncs(update_matrix2::Array{Int,2})
    N = size(update_matrix2)[1]
    update_functions = Array{Function}(undef, N)
    for i in 1:N
        update_functions[i] = (x) -> sign(sum(update_matrix2[:,i].*x))
    end
    return update_functions
end

function checkSSising(update_functions, state)
    N = length(state)
    for i in 1:N
        if state[i] != update_functions[i](state)
            return false
        end
    end
    return true
end

function checkSSnising(update_functions, state)
    N = length(state)
    for i in 1:N
        if state[i] != Int(update_functions[i](state)>0)
            return false
        end
    end
    return true
end

function checkSSnisingTest(update_functions, state)
    N = length(state)
    x = [Int(i>0) for i in state]
    for i in 1:N
        if state[i] != Int(update_functions[i](state)>0)
            x[i] = 0
        else
            x[i] = 1
        end
    end
    return x
end


function asyncIsingNoFunc(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int)
    n_nodes = size(update_matrix,1)
    stateVec = Int[-1,1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_functions = updateFuncs(update_matrix2)
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        flag = 0
        time = 0
        uList = rand(1:n_nodes, nIter)
        for j in 1:nIter
            u = uList[j]
            time = time + 1
            state[u] = update_functions[u](state)
            if iszero(j%5) # check after every five steps,hopefully reduce the time
                if checkSSising(update_functions, state)
                    flag = 1
                    break
                end
            end
        end
        fr = frustration(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
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

function asyncNIsingNoFunc(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int)
    n_nodes = size(update_matrix,1)
    stateVec = Int[0,1]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix2 = 2*update_matrix + Matrix(I, n_nodes, n_nodes)
    update_functions = updateFuncs(update_matrix2)
    @showprogress for i in 1:nInit
        state = rand(stateVec, n_nodes) #pick random state
        init = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
        flag = 0
        time = 0
        uList = rand(1:n_nodes, nIter)
        for j in 1:nIter
            u = uList[j]
            time = time + 1
            state[u] = Int(update_functions[u](state) > 0)
            if iszero(j%5) # check after every five steps,hopefully reduce the time
                if checkSSnising(update_functions, state)
                    flag = 1
                    break
                end
            end
        end
        fr = frustration(state, findnz(sparse(update_matrix)))
        fin = join(["'", join(replace(x -> x == -1 ? 0 : x, state)), "'"])
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