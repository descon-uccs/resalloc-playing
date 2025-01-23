# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:09:58 2025

@author: pbrown2
"""

import scipy.optimize as opt
import numpy as np

# create I
def createISet(n) :
    # this only works for the dual-flavored LP
    I = set()
    for a in range(n+1) :
        for x in range(n+1) :
            for b in range(n+1) :
                if a+x+b>0 and a+x+b<n+1 and (a*x*b == 0 or a+x+b==n):
                    I.add((a,x,b))
    return I

def createCVectorDual(n,w,f,I) :
    # d.v.'s are mu, lambda
    return np.array([1,0])

def createBVectorDual(n,w,f,I) :
    B = np.empty(len(I))
    for i,theta in enumerate(I) :
        a,x,b = theta
        B[i] = -w[b+x]
    return B

def createAMatrixDual(n,w,f,I) :
    A = np.empty([len(I),2])
    for i,theta in enumerate(I) :
        a,x,b = theta
        A[i,0] = -w[a+x]
        A[i,1] = a*f[a+x] - b*f[a+x+1]
    return A
        

if __name__ == "__main__":
    n = 2
    w2 = [0,1,1,1]
    f2 = [0,1,0,0]
    IR2 = createISet(n)
    C = createCVectorDual(n,w2,f2,IR2) 
    B = createBVectorDual(n, w2, f2, IR2)
    A = createAMatrixDual(n, w2, f2, IR2)
    
    result = opt.linprog(C,A,B)
    print(result)