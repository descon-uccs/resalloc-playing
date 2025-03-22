# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 15:17:11 2025

@author: pbrown2
"""


import scipy.optimize as opt
import numpy as np


# create networked I
def createISetNetwork(n) :
    # impractical for n>8 or 9
    I = set()
    for element in range(1,4**n) : # one for each possibility except all-0's
        triples = [None]*n
        for_mod = element
        for j in range(n) :
            config = for_mod % 4
            triple = [0,0,0]
            if config > 0 :
                triple[config-1] = 1
            triples[j] = triple
            for_mod = for_mod//4
        I.add(tuple([tuple(triple) for triple in triples]))
    return I

def checkINetwork(n) :
    # warning: this function can take extremely long to run
    I = createISetNetwork(n)
    sums = {j:0 for j in range(1,n+1)}
    duplicates = 0
    for elem in I :
        accum = 0
        for triple in elem :
            this_sum = sum(list(triple))
            if this_sum > 1 :
                raise Exception('one of the triples has more than a single 1')
            accum += this_sum
        sums[accum] += 1
        for el in I - {elem} :
            if el == elem :
                duplicates += 1
    return sums, duplicates
        
def AMatrix(n,w,f,I) :
    # f is length-n
    A = np.empty([n+1,len(I)])
    for i,theta in enumerate(I) :
        for j in range(n) :
            aj = theta[j][0]
            bj = theta[j][2]
            Atj = sum([theta[k][0] + theta[k][1] for k in range(n)]) # note: with non-complete network, this needs to be changed!
            A[j,i] = -(aj*f[j][Atj] - bj*f[j][Atj+1])
        A[n,:] = A[:n,:].sum(0)
    return A

def Bvector(n,w,f,I,S=0,mode='unconstrained') :
    '''
    Modes:
        - 'unconstrained': stability can be allocated to agents in any way
        - 'one': all stability is allocated to a single agent
        - 'even': stability is split evenly among the agents
    '''
    if mode=='unconstrained' :
        B = np.array([[0]*n + [-S]]).T
    if mode == 'one' :
        B = np.array([[-S]+[0]*(n-1) + [-S]]).T
    if mode == 'even' :
        B = np.array([[-S/n]*n + [-S]]).T
    return B

def AMatrixEq(n,w,f,I) :
    A = np.empty([1,len(I)])
    for i,theta in enumerate(I) :
        At = sum([theta[k][0] + theta[k][1] for k in range(n)])
        A[0,i] = w[At]
    return A

def BvectorEq(n,w,f,I) :
    return np.array([1])

def Cvector(n,w,f,I) :
    C = np.empty(len(I))
    for i,theta in enumerate(I) :
        Bt = sum([theta[k][2] + theta[k][1] for k in range(n)])
        C[i] = -w[Bt]
    return C

        
if __name__== "__main__" :
    n = 7
    I = createISetNetwork(n)
    w = [0] + [1 for _ in range(n)]
    fes = [0] + [w[i]/i for i in range(1,len(w))]+[0]
    fmc = [0] + [w[i] - w [i-1] for i in range(1,len(w))] + [0]
    # fmc = [0,1] + [0 for i in range(1,n+1)]
    
    
    S=0.1
    mode='even'
    
    f = [fes]*(n)
    A = AMatrix(n,w,f,I)
    B = Bvector(n, w, f, I,S=S,mode=mode)
    Aeq = AMatrixEq(n,w,f,I)
    Beq = BvectorEq(n,w,f,I)
    C = Cvector(n, w, f, I)
    
    result = opt.linprog(C,
                         A_ub=A,
                         b_ub=B,
                         A_eq=Aeq,
                         b_eq=Beq)
    print(-1/result.fun)
    
    