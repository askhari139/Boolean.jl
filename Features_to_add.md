- Code configuration files to deal with the plaethora of options. Add also a tree to demonstrate the simulation order. Maybe include a network object that holds the parameters?
- New methods
    1. nearest level rather than fixed intervals for levels. This might basically bias the system largely towards the extremes (since they have a larger range : from 1 to inf or -inf to -1.). But, that is if all values in the number line have equal chance of getting picked. 
    2. Multilevel nIsing formalism : 0, 0.5, 1 instead of -1, -0.5, 0.5, 1
    3. Multilevel as -n to n instead of -1 to 1
    4. STG driven simulations
    5. Synchronous updates
    6. add a threshold parameter in nIsing, default being 0
- Partial turn-off and partial multilevel
- State transition graph
- Basin parameters
    1. list of hamming distances of init conditions
- Initial conditions that converge to multiple steady states
- Multi-level formalism: update the nodes only if there is a change.
- Implement QSSA in boolean formalism
- Implement logic based functions
- Trajectories for specific states, preferably pick the states from the simulation outputs
- Derrida function calculation
- Logical rule modelling:
    - Using prime implicants and trap spaces to accelerate the computation
    
