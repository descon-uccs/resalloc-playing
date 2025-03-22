# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 15:17:11 2025

@author: pbrown2
"""


import scipy.optimize as opt
import numpy as np

# create I
def createISetPrimal(n) :
    I = set()
    for a in range(n+1) :
        for x in range(n+1) :
            for b in range(n+1) :
                if a+x+b>0 and a+x+b<n+1 :
                    I.add((a,x,b))
    return I

def createCVectorPrimal(n,w,f,I) :
    C = np.empty(len(I))
    for i,theta in enumerate(I) :
        a,x,b = theta
        C[i] = -w[b+x]
    return C

def createBVectorPrimal(n,w,f,I,S=0) :
    return np.array([-S,1])

def createAMatrixPrimal(n,w,f,I) :
    A = np.empty([2,len(I)])
    for i,theta in enumerate(I) :
        a,x,b = theta
        A[0,i] = -(a*f[a+x] - b*f[a+x+1])
        A[1,i] = w[a+x]
    return A
        

if __name__ == "__main__":
    n = 4
    w = [0] + [1 for _ in range(n)]
    fes = [0] + [w[i]/i for i in range(1,len(w))]+[0]
    fmc = [0] + [w[i] - w [i-1] for i in range(1,len(w))] + [0]
    # fmc = [0,1] + [0 for i in range(1,n+1)]
    I = createISetPrimal(n)
    
    f = fes
    C = createCVectorPrimal(n,w,f,I)
    B = createBVectorPrimal(n, w, f, I,S=0)
    A = createAMatrixPrimal(n, w, f, I)
    
    result = opt.linprog(C,A_eq=A,b_eq=B) # doing it all with equality constraints to hit S
    print(-1/result.fun)