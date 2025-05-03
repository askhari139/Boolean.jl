import numpy as np
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


def f_from_expression(expr):
    # split the expression by "=" and remove leading and trailing spaces
    expr = expr.strip()
    expr = expr.replace('AND', '&').replace('NOT', '!').replace('OR', '|').replace('!','~ ')
    expr = expr.replace('(',' ( ').replace(')',' ) ')
    expr_split = expr.split(' ')
    
    var = []
    dict_var = dict()
    n_var = 0
    for i,el in enumerate(expr_split):
        if el not in ['',' ','(',')','and','or','not','AND','OR','NOT','&','|','~','+','-','*','%','>','>=','==','<=','<'] and not el.isdigit():
            try:
                new_var = dict_var[el]
            except KeyError:
                new_var = 'x[%i]' % n_var
                dict_var.update({el:new_var})
                var.append(el)
                n_var += 1
            expr_split[i] = new_var
    expr = ' '.join(expr_split)
    F = []
    for x in itertools.product([0, 1], repeat = n_var):
        F.append(int(eval(expr))) #x is used here "implicitly"
    return F,var

def is_monotonic(F,GET_DETAILS=False, non_monotonic = False):
    n=int(np.log2(len(F)))
    F = np.array(F)
    monotonic = []
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        diff = F[dummy==1]-F[dummy==0]
        min_diff = min(diff)
        max_diff = max(diff)
        if min_diff==0 and max_diff==0:
            monotonic.append('0')
        elif min_diff==-1 and max_diff==1:
            monotonic.append('0')
        elif min_diff>=0 and max_diff==1:
            monotonic.append('1')            
        elif min_diff==-1 and max_diff<=0:
            monotonic.append('2')   
    if GET_DETAILS:
        return ('0' not in monotonic,monotonic)
    else:
        return '0' not in monotonic

def f_from_file(booleanRules, n_max = 100):
    fl = []
    with open(booleanRules) as f_in:
        for line in f_in:
            if line[0] != '#':
                fl.append(line.strip())
    if (len(fl) > n_max):
        return 0
    table = dict()
    F, I, degree, var, constants = text_to_BN(booleanRules, max_N = n_max)
    variables = var + constants
    varList = []
    for inputs in I:
        varList.append([variables[i] for i in inputs])
    for i in range(len(var)):
        table.update({var[i]:[F[i], I[i]]})
    # for line in fl:
    #     line = line.split('=')
    #     f,variables = f_from_expression(line[1])
    #     table.update({line[0].strip():[f, variables]})
    return table

def boolean_to_topo(BooleanRules, randomize=False, effectivity=False, n_max = 100):
    TruthTable = f_from_file(BooleanRules, n_max = n_max)
    if TruthTable == 0:
        print('File too long')
        return 0
    topoFile = []
    output_file = BooleanRules.split('.')[0]
    topoFile.append('Source Target Type')
    # For each value in TruthTable, find monotonicity
    for key in TruthTable:
        f = TruthTable[key][0]
        k = int(np.log2(len(f)))
        if (randomize):
            random.shuffle(f)
        MONOTONIC, regulation_types = is_monotonic(f, GET_DETAILS=True)
        
        for i in range(len(TruthTable[key][1])):
            if regulation_types[i] == "0":
                continue
            node = TruthTable[key][1][i]
            topoFile.append(node + ' ' + key + ' ' + regulation_types[i])
        
    with open(output_file + '.topo', 'w') as f:
        for line in topoFile:
            f.write(line + '\n')
    return 1

def bools_to_topo():
    # list of all txt files
    import os
    files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.txt')]
    files.sort()
    # loop through each file
    for file in files:
        print(file)
        r = boolean_to_topo(file, n_max = 100)
        if (r == 0):
            print('File too long')
            continue
        # read each line of the topo file
        with open(file.split('.')[0] + '.topo', 'r') as f:
            lines = f.readlines()
        # loop through each line and check if it ends with a number
        linesUpd = [lines[0]]
        for line in lines:
            if line.endswith('1\n') or line.endswith('2\n'):
                # if it does, remove the last character
                linesUpd.append(line)
        # write the updated lines to a new file
        with open(file.split('.')[0] + '.topo', 'w') as f:
            for line in linesUpd:
                f.write(line)


def F_to_topo(F, output_file, randomize=False):
    TruthTable = F
    topoFile = []
    topoFile.append('Source Target Type')
    # effectivity = []
    # effectivity.append('Source,Target,Effectiveness')
    # For each value in TruthTable, find monotonicity
    for key in TruthTable:
        f = TruthTable[key][0]
        k = int(np.log2(len(f)))
        if (randomize):
            random.shuffle(f)
        MONOTONIC, regulation_types = is_monotonic(f, GET_DETAILS=True)
        # f_cana = cana.boolean_node.BooleanNode(inputs=range(k),k=k,outputs=f)
        # eff = f_cana.edge_effectiveness()
        for i in range(len(TruthTable[key][1])):
            if regulation_types[i] == "0":
                continue
            node = TruthTable[key][1][i]
            topoFile.append(node + ' ' + key + ' ' + regulation_types[i])
            # effectivity.append(node + ',' + key + ',' + str(eff[i]))

    #write topoFile to output_file
    with open(output_file + '.topo', 'w') as f:
        for line in topoFile:
            f.write(line + '\n')

    # #write effectivity to output_file
    # with open(output_file + '.eff', 'w') as f:
    #     for line in effectivity:
    #         f.write(line + '\n')
    return
