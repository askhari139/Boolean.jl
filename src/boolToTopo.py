#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:17:59 2022

@author: ckadelka
"""

import numpy as np
import itertools
import random
import cana
import cana.boolean_node

## Functions to determine the type of regulation of a boolean function
# Input: Boolean function in the form of a list of 0s and 1s
# Output: Boolean function in the form of a list of 0s and 1s
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
    for line in fl:
        line = line.split('=')
        f,variables = f_from_expression(line[1])
        table.update({line[0].strip():[f, variables]})
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
