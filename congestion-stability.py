# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 15:57:23 2025

@author: pbrown2
"""



import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt

# create I
def createISetPrimal(n,m) :
    I = set()
    for x in range(n+1) :
        for y in range(n+1) :
            for z in range(n+1) :
                for j in range(m) :
                    if x+y+z>0 and x+y+z<n+1 :
                        I.add((x,y,z,j))
    return I

def createCVectorCongestion(n,b,f,I) :
    # b is a list of basis functions or lambdas
    C = np.empty(len(I))
    for i,theta in enumerate(I) :
        x,y,z,j = theta
        C[i] = b[j](x+z)*(x+z)
    return C

def createBVectorCongestion(S=0) :
    return np.array([-S])

def createBeqVectorCongestion() :
    return np.array([1])

def createAMatrixCongestion(n,b,f,I) :
    # Nash constraint
    A = np.empty([1,len(I)]) 
    for i,theta in enumerate(I) :
        x,y,z,j = theta
        A[0,i] = f[j](x+y)*y - f[j](x+y+1)*z
    return A

def createAeqMatrixCongestion(n,b,f,I) :
    # Objective normalization constraint
    A = np.empty([1,len(I)])
    for i,theta in enumerate(I) :
        x,y,z,j = theta
        A[0,i] = b[j](x+y)*(x+y)
    return A
        
def getPoA(n,b,f,S) :
    
    I = createISetPrimal(n,len(b))
    C = createCVectorCongestion(n,b,f,I)
    B = createBVectorCongestion(S)
    B_eq = createBeqVectorCongestion()
    A = createAMatrixCongestion(n,b,f,I)
    A_eq = createAeqMatrixCongestion(n,b,f,I)
    
    result = opt.linprog(C,A,B,A_eq,B_eq,method='simplex')
    return 1/result.fun

def PoASmoothness(lamb,mu,S) :
    return lamb/(1-mu+S)


def findLambdaMu(d) :
    # this is a grid search implementation of Lemma 5.8 in "Exact PoA for Polynomial Congestion Games" by Aland et al.
    
    expr = lambda x,mu: pow((x+1),d)-mu*pow(x,d+1)
    
    minwrtmu = np.inf
    minmu = 0.01
    for mu in np.linspace(.01,.99,99) :
        maxwrti = -np.inf
        for x in range(1000) :
            maxwrti = max(maxwrti, expr(x,mu)/(1-mu))
        if maxwrti < minwrtmu :
            minwrtmu = maxwrti
            minmu = mu
    return minwrtmu*(1-minmu), minmu


if __name__ == "__main__":
    n = 4
    
    b = [lambda x: 1, lambda x: x, lambda x: x*x]
    
    # ftest = [lambda x: 1, lambda x: x+1]
    fself = b
    fmc = [lambda x: b[0](x)+(x-1)*(b[0](x)-b[0](x-1)), lambda x: b[1](x)+(x-1)*(b[1](x)-b[1](x-1))]
    
    f=fself
    SS=np.linspace(0,1)
    
    d = len(b)-1
    lamb,mu = findLambdaMu(d)
    PoAs = [getPoA(n,b,f,S) for S in SS]
    PoAsLambdaMu = [PoASmoothness(lamb,mu,S) for S in SS]
    plt.plot(SS,PoAs)
    plt.plot(SS,PoAsLambdaMu)