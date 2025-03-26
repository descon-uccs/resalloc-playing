# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 15:17:11 2025

@author: pbrown2
"""


import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from math import factorial
import os

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

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

class PoAPlotter :
    def __init__(self,n,wtype='setcover',ftype='mc') :
        self.n = n
        self.I = createISetNetwork(n)
        if wtype=='setcover':
            self.w = [0] + [1 for _ in range(n)]
        w = self.w
        if ftype == 'mc':
            fmc = [0] + [w[i] - w[i-1] for i in range(1,len(w))] + [0]
            self.f = [fmc]*n
        elif ftype == 'es':
            fes = [0] + [w[i]/i for i in range(1,len(w))]+[0]
            self.f = [fes]*n
        elif ftype == 'opt':
            thesum = lambda j : sum([1/factorial(i) for i in range(0,j)])
            fopt = [0] + [factorial(j-1)/(np.e-1)*(np.e-thesum(j)) for j in range(1,n+3)] # longer than we need, it's ok
            self.f = [fopt]*n
        self.smallfontsize = 12
        self.medfontsize = 14
        self.largefontsize = 16
        self.linewidth = 3
    
    def PoA(self,S=0,mode='unconstrained') :
        I = self.I
        n = self.n
        w = self.w
        f = self.f
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
        return -1/result.fun
    
    def plotPoA(self,Srange,mode='unconstrained',fignum=13,label=None,ax=None) :
        # Srange is a tuple defining a closed interval for stability margin
        SS = list(np.linspace(Srange[0], Srange[1], 50))
        PoAs = []
        for S in SS :
            PoAs.append(self.PoA(S=S,mode=mode))
        
        if not ax :
            fig = plt.figure(num=fignum)
            fig.clf()
            ax = fig.add_subplot(1,1,1)
        ax.plot(SS,PoAs,label=label, linewidth=self.linewidth)
        if label :
            ax.legend()
        ax.set_xlim(left=0)
        ax.set_ylim([0.5,1])
    
        
        # Styling
        ax.set_title("Price of Anarchy vs. S", fontsize=self.largefontsize)
        ax.set_xlabel("S", fontsize=self.medfontsize)
        ax.set_ylabel("Price of Anarchy", fontsize=self.medfontsize)
        
        # Gridlines
        ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        
        # Legend
        if label:
            ax.legend(fontsize=self.smallfontsize)
        
        # Optional: nicer tick label formatting
        ax.tick_params(axis='both', which='major', labelsize=self.smallfontsize)
        
        return ax
    
    
def save_ax_as_pdf(ax, filename="my_plot.pdf"):
    directory = 'figures'
    fig = ax.figure
    os.makedirs(directory, exist_ok=True)
    fig.savefig(directory+'/'+filename, format="pdf", bbox_inches="tight")


def createPlots(save=False) :
    # call this function to create all the plots we need (and more!)
    
    # plots with unconstrained stability margin, ES (various n) and MC (n=2)
    fignum = 35
    plotter2MC = PoAPlotter(2,ftype='mc')
    ESplotters = [PoAPlotter(n,ftype='es') for n in range(2,7)]
    ESax = None
    for plotter in ESplotters :
        n = plotter.n
        label = 'ES, n=' + str(n)
        if ESax :
            plotter.plotPoA([0,(n-1)/n],mode='unconstrained',ax=ESax,label=label)
        else :
            ESax = plotter.plotPoA([0,(n-1)/n],mode='unconstrained',fignum=fignum,label=label)
    ESax.set_xlim([0,1])
    plotter2MC.plotPoA([0,1],ax=ESax,label='MC, n=2')
    if save :
        save_ax_as_pdf(ESax,'Figure_1.pdf')
    
    # conjecture: the ES PoA is given by 1/((2n-1)/n - S). Matches nicely what we have here!
    
    # can uncomment the following to get a "Gairing" plot, but it's not very interesting. Actually it sits *almost* exactly between the n=2 and n=3 ES traces.
    # plotter2Opt = PoAPlotter(2,ftype='opt')
    # plotter2Opt.plotPoA([0,.6],ax = ESax,label = 'OPT, n=2')
    
    
    # plots for ES, individual values of n, comparing unconstrained SM with putting all the SM on a single player
    axOne = None
    fignum = 36
    for plotter in ESplotters :
        n = plotter.n
        label = 'ES, n='+str(n)
        uncon = ', unconstrained'
        one = ', one'
        axOne = plotter.plotPoA([0,(n-1)/n],mode='unconstrained',fignum=fignum+n,label=label+uncon)
        plotter.plotPoA([0,0.6472],mode='one',ax=axOne,label=label+one)
        if save :
            save_ax_as_pdf(axOne,str(fignum+n)+'.pdf')
        
        
if __name__== "__main__" :
    
    
    createPlots(save=True)
    
    
    