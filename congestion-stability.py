# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 15:57:23 2025

@author: pbrown2
"""



import scipy.optimize as opt
import numpy as np

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
    n = 2
    
    b = [lambda x: 1, lambda x: x]
    # fmc = [0,1] + [0 for i in range(1,n+1)]
    I = createISetPrimal(n,2)
    
    f = fes
    C = createCVectorPrimal(n,w,f,I)
    B = createBVectorPrimal(n, w, f, I,S=0)
    A = createAMatrixPrimal(n, w, f, I)
    
    result = opt.linprog(C,A,B)
    print(-1/result.fun)