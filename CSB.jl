function csb(update_matrix::Array{Int,2},
    nInit::Int, nIter::Int, timeStep::Float64)
    n_nodes = size(update_matrix,1)
    ls = collect(1:nLevels)
    stateVec = [-1*reverse(ls)/nLevels; ls/nLevels]
    initVec = []
    finVec = []
    flagVec = []
    frustVec = []
    timeVec = []
    # states_df = DataFrame(init = String[], fin = String[], flag = Int[])
    update_matrix = float(update_matrix)
    for i in 1:nInit
        state = rand(stateVec, n_nodes)
        state = stateConvert(state, nLevels = nLevels)
        push!(initVec, state)
        state = state[2:end-1]
        state = [parse(Int, x) for x in split(state, ",")]
        state = stateConvert(state, nLevels = nLevels)
        push!(finVec, state)
        flag = 0
        for j in 1:nIter
            state = update_matrix*state
            state = stateConvert(state, nLevels = nLevels)
            if state in initVec
                flag = 1
                break
            end
        end
        push!(flagVec, flag)
        push!(frustVec, shubhamFrust(state, nonZeros))
        push!(timeVec, j)
    end
    states_df = DataFrame(init = initVec, fin = finVec, flag = flagVec, 
        frust = frustVec, time = timeVec)
    return states_df
end