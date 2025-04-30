import numpy as np
import pandas as pd
import canalizing_function_toolbox_v16 as can
import load_database13 as db
import random

import itertools
import json

def find_all_indices(arr,el):
    '''
    Given a list arr, this function returns a list of all the indices i where arr[i]==el.
    If el not in arr, it raises a ValueError.
    '''
    res=[]
    for i,a in enumerate(arr):
        if a==el:
            res.append(i)
    if res==[]:
        raise ValueError('The element is not in the array at all')
    return res

def text_to_BN(textfile,separator_var_func="=",original_not="NOT",original_and="AND",original_or="OR",new_not=" not ",new_and=" and ",new_or=" or ",max_degree=15,TREATMENT_OF_CONSTANTS=1,max_N=10000):
    '''
    This function takes as input a textfile in directory folder, 
    where each line describes the regulation of one gene, specified to the left, 
    separated by 'separator_var_func' from its update rules to the right.
    
    The function outputs a Boolean network model, specified by
        F: a list of N+C lists of length 2^(n_i) 
        where N is the number of genes, C the number of constants (i.e., external parameters) 
        and n_i the number of regulators per gene (n_i=1 for all constants by definition).
        F describes the update rules of the BN, stored as the right-hand side of a truth table
        
        I: a list of N lists of length n_i where N is the number of genes, n_i the number of regulators per gene. 
        I describes the regulators of each gene (as indices 0, 1, ..., n-1), i.e., the wiring diagram of the BN (the adjacency matrix of a directed graph)
        
        degree = [n_1, n_2, ..., n_N]: a list of length N, describing the number of regulators of each gene
        
        var: a list of the names of the N genes
        
        constants: a list of the names of the C constants or external parameters
    
    Inputs:
        original_not, original_and, original_or: Boolean logic operators used in the updated rules of the provided text file
        
        new_not, new_and, new_or: Boolean logic operators to be used (default: Python interpretable Boolean logic operator)
        
        max_degree: (default 15) update rules with more than max_degree regulators are not computed, instead an empty list [] is added as a placeholder
        
        max_N: (default 10,000) text files with more than max_N rows yield an AssertionError
            
        TREATMENT_OF_CONSTANTS: Ternary choice:
            0: constants (i.e., external parameters) are not added to the BN, yields a BN that cannot be dynamically evaluated and causes errors unless only the degree distribution and update rules are studied,
            1: (default) constants are added as self-regulatory nodes into the network, which is then by definition not strongly-connected,
            2: (not implemented yet) multiple models are returned, one for each combination of constants, the constants are not included as nodes but instead the update rules are simplified

    Example of an input file:
        A = NOT B
        B = A OR C
        C = E OR (A AND (NOT B))
    
    Output with TREATMENT_OF_CONSTANTS==1 (default):
        F = [[1,0],
             [0,1,1,1],
             [0,1,0,1,1,1,0,1],
             [0,1]]
        I = [[1],
             [0,2],
             [0,1,3],
             [3]]
        degree = [1,2,3,1]
        var = ['A','B','C']
        constants = ['E']
        
    Output with TREATMENT_OF_CONSTANTS==0:
        F = [[1,0],
             [0,1,1,1],
             [0,1,0,1,1,1,0,1]]
        I = [[1],
             [0,2],
             [0,1,3]]
        degree = [1,2,3]
        var = ['A','B','C']
        constants = ['E']    
    '''
    
    f = open(textfile,'r')
    text = f.read()
    text = text.replace('\t',' ').replace('(',' ( ').replace(')',' ) ')
    tvec = text.splitlines()
    f.close()
    
    #Deleting empty lines
    while '' in tvec:
        tvec.remove('')
        
    n=len(tvec)
    assert n<=max_N,'n='+str(n)+' > max_N='+str(max_N)
    var=["" for i in range(n)]
    
    #Determining Variables, var contains the order of variables (may differ from original order)
    for i in range(n):
        var[i]=tvec[i][0:tvec[i].find(separator_var_func)].replace(' ','')

    constants_and_variables = []
    for line in tvec:
        linesplit = line.split(' ')
        for el in linesplit:
            if el not in ['(',')','+','*','1',separator_var_func,original_not,original_and,original_or,'',' ']:
                constants_and_variables.append(el)
    constants = list(set(constants_and_variables)-set(var))
        
    dict_variables_and_constants = dict({original_not:new_not,original_and:new_and,original_or:new_or})
    dict_variables_and_constants.update(dict(list(zip(var,["x[%i]" % i for i in range(len(var))]))))
    dict_variables_and_constants.update(list(zip(constants,["x[%i]" % i for i in range(len(var),len(set(constants_and_variables)))]))) #constants at end

    for i,line in enumerate(tvec):
        linesplit = line.split(' ')
        for ii,el in enumerate(linesplit):
            if el not in ['(',')','+','*','1',separator_var_func,new_not.strip(' '),new_and.strip(' '),new_or.strip(' '), '',' ']:
                linesplit[ii] = dict_variables_and_constants[el]
        tvec[i] = ' '.join(linesplit)
    #       
    for ind in range(n):
        tvec[ind]=tvec[ind][tvec[ind].find(separator_var_func)+len(separator_var_func):]
        #tvec[ind]="("+tvec[ind]+")"

    #create I, the list of essential variables
    I = []
    tvec_mod = []
    for i in range(n):
        indices_open = find_all_indices(tvec[i],'[')
        indices_end = find_all_indices(tvec[i],']')
        dummy = np.sort(np.array(list(map(int,list(set([tvec[i][(begin+1):end] for begin,end in zip(indices_open,indices_end)]))))))
        I.append( dummy )
        dict_dummy = dict(list(zip(dummy,list(range(len(dummy))))))
        tvec_dummy = tvec[i][:]
        for el in dummy:
            tvec_dummy = tvec_dummy.replace('[%i]' % el,'[%i]' % dict_dummy[el]) #needs to be done with an ascending order in dummy
        tvec_mod.append(tvec_dummy)
    
    degree = list(map(len,I))
        
    F = []
    for i in range(n):
        f = np.array([],dtype=int)
        if degree[i]<=max_degree:
            X = list(itertools.product([0, 1], repeat = degree[i]))
            for j in range(2**degree[i]):
                x = X[j] #x is used "implicitly" in the next line!!
                f = np.append(f,eval(tvec_mod[i])%2) #x is used here "implicitly"
        F.append(f)
        
    
    if TREATMENT_OF_CONSTANTS==1:
        for i in range(len(constants)):
            F.append(np.array([0,1]))
            I.append(np.array([len(var)+i]))
            degree.append(1)
    assert TREATMENT_OF_CONSTANTS in [0,1],'TREATMENT_OF_CONSTANTS must be 0 or 1 (default). TREATMENT_OF_CONSTANTS==2, yielding 2^C models for each combination of inputs to the C constants is not yet implemented.'
    result = {
        "F": [f.tolist() for f in F],
        "I": [i.tolist() for i in I],
        "degree": degree,
        "variables": var,
        "constants": constants
    }
    outputfile = textfile.replace('.txt','.json')
    with open(outputfile, 'w') as out:
        json.dump(result, out, indent=2)
    return F, I, degree, var, constants

def load_database(folders = ['update_rules_cell_collective/'],separator_var_func="=",original_not="NOT",original_and="AND",original_or="OR",new_not=" not ",new_and=" and ",new_or=" or ",max_degree=15,max_N=10000):
    '''
    This function takes as input a list of directories, in which it searches for Boolean networks stored as tabular models or in text format.
    
    The function outputs a list of all M Boolean network models it was able to convert to a standardized format.
    In this format, each BN is specified by
        F: a list of N+C lists of length 2^(n_i) 
        where N is the number of genes, C the number of constants (i.e., external parameters) 
        and n_i the number of regulators per gene (n_i=1 for all constants by definition).
        F describes the update rules of the BN, stored as the right-hand side of a truth table
        
        I: a list of N lists of length n_i where N is the number of genes, n_i the number of regulators per gene. 
        I describes the regulators of each gene (as indices 0, 1, ..., n-1), i.e., the wiring diagram of the BN (the adjacency matrix of a directed graph)
        
        degree = [n_1, n_2, ..., n_N]: a list of length N, describing the number of regulators of each gene
        
        degree_essential = [e_1, e_2, ..., e_N]: a list of length N, describing the number of ESSENTIAL regulators of each gene, 
        where a regulator is essential if has an impact on the function
        
        var: a list of the names of the N genes
        
        constants: a list of the names of the C constants or external parameters
        
    The actual outputs of the function are:
        Fs: a list of M F's
        
        Is: a list of M I's
        
        degrees: a list of M degree's
        
        degrees_essential: a list of M degree_essential's
            
        variabless: a list of M var's
        
        constantss: a list of M constants's 
        
        models_loaded: a list of M filenames of the models that were successfully converted
        
        models_not_loaded: a list of filenames of the models that were NOT successfully converted
    
    Optional inputs:
        original_not, original_and, original_or: Boolean logic operators used in the updated rules of the provided text file
        
        new_not, new_and, new_or: Boolean logic operators to be used (default: Python interpretable Boolean logic operator)
        
        max_degree: (default 15) update rules with more than max_degree regulators are not computed, instead an empty list [] is added as a placeholder
        
        max_N: (default 10,000) text files with more than max_N rows yield an AssertionError
            
        TREATMENT_OF_CONSTANTS: Ternary choice:
            0: constants (i.e., external parameters) are not added to the BN, yields a BN that cannot be dynamically evaluated and causes errors unless only the degree distribution and update rules are studied,
            1: (default) constants are added as self-regulatory nodes into the network, which is then by definition not strongly-connected,
            2: multiple models are returned, one for each combination of constants, the constants are not included as nodes but instead the update rules are simplified
    '''

    Fs,Is,degrees,variabless,constantss,degrees_essential = [],[],[],[],[],[]
    models_loaded,models_not_loaded = [],[]
    
    for folder in folders:
        for fname in os.listdir(folder):
            if fname.endswith('tabular.txt'): #first check if it is a model that is already available in tabular form
                try:
                    textfile = fname
                    F,I,degree,variables, constants = load_tabular_model(folder,fname,max_N=max_N)
                    print(textfile,'converted')
                except:
                    models_not_loaded.append(textfile)
                    print()
                    print(textfile,'failed')
                    print()
                    continue
            elif fname.endswith('.txt'):
                try:
                    textfile = fname
                    F,I,degree,variables, constants = text_to_BN(folder,textfile,max_degree=max_degree,max_N=max_N)
                    print(textfile,'converted')
                except:
                    models_not_loaded.append(textfile)
                    print()
                    print(textfile,'failed')
                    print()
                    continue
            else:
                continue
            #if len(constants)>0:
            #    print(textfile,'has %i constants' % len(constants))
                    
            models_loaded.append(textfile)            
            Fs.append(F)
            Is.append(I)
            degrees.append(degree)
            degrees_essential.append([can.get_number_essential_variables(f) for f in F])
            variabless.append(variables)
            constantss.append(constants)
    return [Fs,Is,degrees,degrees_essential,variabless,constantss,models_loaded,models_not_loaded]


def bin2dec(binary_vector):
    """
    Convert a binary vector to an integer.

    Parameters:
    binary_vector (list): List containing binary digits (0 or 1).

    Returns:
    int: Integer value converted from the binary vector.
    """
    binary_string = ''.join(str(bit) for bit in binary_vector)
    return int(binary_string, 2)


def dec2bin(integer_value, num_bits):
    """
    Convert an integer to a binary vector.

    Parameters:
    integer_value (int): Integer value to be converted.
    num_bits (int): Number of bits in the binary representation.

    Returns:
    list: List containing binary digits (0 or 1).
    """
    binary_string = bin(integer_value)[2:].zfill(num_bits)
    return [int(bit) for bit in binary_string]


def update(F, I, N, X):
    Fx = np.zeros(N, dtype = int)
    for i in range(N):
        # print(i)
        Fx[i] = F[i][bin2dec(X[I[i]])]
    return Fx

def update_single_node(f,states_regulators):
    return f[bin2dec(states_regulators)]

def num_of_steady_states_asynchronous(F,I,N,nsim=500,EXACT = False, left_side_of_truth_table = [], initial_sample_points = [],search_depth=1000,SEED=-1,DEBUG=False): 
    if EXACT and left_side_of_truth_table == []:
        left_side_of_truth_table = list(map(np.array,list(itertools.product([0, 1], repeat = N))))

    sampled_points = []
    
    assert initial_sample_points==[] or not EXACT,"sample points povided but with option EXACT the entire state space is computed (initial sample points ignored)"
    
    if SEED == -1:
        SEED = int(random.random()*2**31)
    
    np.random.seed(SEED)
    
    dictF = dict()
    steady_states = []
    basin_sizes = []
    steady_state_dict = dict()   
    for iteration in range(nsim if not EXACT else 2**N):
        if EXACT:
            x = left_side_of_truth_table[iteration]
            xbin = iteration
        else:
            if initial_sample_points==[]: #generate random initial states on the fly
                x = np.random.randint(2, size = N)
                xbin = bin2dec(x)
                sampled_points.append(xbin)
            else:                
                x = initial_sample_points[iteration]
                xbin = bin2dec(x)
        
        if DEBUG:
            print(iteration,-1,-1,False,xbin,x)        
        for jj in range(search_depth): #as long as we haven't reached an attractor state, keep updating
            FOUND_NEW_STATE = False
            try: # check if this state is a known steady state
                index_ss = steady_state_dict[xbin]
            except KeyError: #update the state asynchrnonously until a steady state is reached
                update_order_to_try = np.random.permutation(N)
                for i in update_order_to_try:
                    try:
                        fxbin = dictF[(xbin,i)]
                        if fxbin!=xbin:
                            FOUND_NEW_STATE = True
                            x[i] = 1-x[i]
                    except KeyError:
                        fx_i = update_single_node(F[i],x[I[i]])
                        if fx_i > x[i]:
                            fxbin = xbin + 2**(N-1-i)
                            x[i] = 1
                            FOUND_NEW_STATE = True
                        elif fx_i < x[i]:
                            fxbin = xbin - 2**(N-1-i)
                            x[i] = 0
                            FOUND_NEW_STATE = True
                        else:
                            fxbin = xbin
                        dictF.update({(xbin,i):fxbin})
                    if FOUND_NEW_STATE:
                        xbin = fxbin
                        break
                if DEBUG:
                    print(iteration,jj,i,FOUND_NEW_STATE,xbin,x)
            if FOUND_NEW_STATE == False: #found a steady state in the asynchronous state space
                try: # check if this state is a known steady state
                    index_ss = steady_state_dict[xbin] #returns a KeyError if we haven't identified fxbin as a staedy before
                    basin_sizes[index_ss] += 1
                    break
                except KeyError:
                    steady_state_dict.update({xbin:len(steady_states)})
                    steady_states.append(xbin)
                    basin_sizes.append(1)
                    break
        if DEBUG:
            print()
    if sum(basin_sizes)<(nsim if not EXACT else 2**N):
        print('Warning: only %i of the %i tested initial conditions eventually reached a steady state. Try increasing the search depth. It may however also be the case that your asynchronous state space contains a limit cycle.' % (sum(basin_sizes),nsim if not EXACT else 2**N))
    return (steady_states,len(steady_states),basin_sizes,steady_state_dict,dictF,SEED,initial_sample_points if initial_sample_points!=[] else sampled_points)

