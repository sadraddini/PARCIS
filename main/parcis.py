#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 14:07:39 2018

@author: sadra
"""

"""
    Here I try a different method
"""

from gurobipy import *
import numpy as np

class system:
    def __init__(self,n,m):
        self.n=n
        self.m=m
        self.A=np.empty((n,n))
        self.B=np.empty((n,m))
        self.c=np.empty((n,1))
        self.x=np.empty((n,1))
        self.u=np.empty((m,1))
        self.Pi={} # P matrix, n_pi*n
        self.H={} #
        self.h={} #
        self.F={}
        self.f={}
        self.W={} # Matrix characterizing W
        self.phi={}
        self.theta={}
        self.x_bar=np.empty((n,1))
        self.u_bar=np.empty((m,1))

    def RCI(self,K=3):
        model=Model("RCI_set")
        (n,n_w)=self.W.shape
        if n!=self.n:
            raise(ValueError("the number of rows in W=%d does not match n=%d"%(n,self.n)))
        T=K*n_w
        self.Pi=pi_cube(T) # This is the box assumpion
        phi=np.empty((self.n,T),dtype='object')
        theta=np.empty((self.m,T),dtype='object')
        x_bar=np.empty((self.n,1),dtype='object')
        u_bar=np.empty((self.m,1),dtype='object')
        for row in range(self.n):
            x_bar[row]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
            for column in range(T):
                phi[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for row in range(self.m):
            u_bar[row]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
            for column in range(T):
                theta[row,column]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        model.update()
        # dynamics: A*phi+B*theta: first n_w columns become zero, rest with W become H!
        for row in range(self.n):
            for column in range(T):
                exp=LinExpr()
                for k in range(self.n):
                    exp.add(self.A[row,k]*phi[k,column])
                for k in range(self.m):
                    exp.add(self.B[row,k]*theta[k,column])
                if column<n_w:
                    model.addConstr(exp==0)
                else:
                    model.addConstr(exp==phi[row,column-n_w])
        for row in range(self.n):
            for column in range(n_w):
                model.addConstr(phi[row,T-n_w+column]==self.W[row,column])
        # Fixed point: A*x_bar+B*u_bar+c=x_bar
        for row in range(self.n):
            exp=LinExpr()
            for k in range(self.n):
                exp.add(self.A[row,k]*x_bar[k,0])
            for k in range(self.m):
                exp.add(self.B[row,k]*u_bar[k,0])
            model.addConstr(exp+self.c[row,0]==x_bar[k,0])
        # Everything inside H:
        subset(model,phi,self.Pi,self.H,self.h,x_bar)
        subset(model,theta,self.Pi,self.F,self.f,u_bar)
        # Solve the optimization problem
        model.optimize()
        self.phi=valuation(phi)
        self.theta=valuation(theta)
        self.x_bar=valuation(x_bar)
        self.u_bar=valuation(u_bar)


def subset(model,G,Pi,H,h,x):
    """
    Description: Add Farkas lemma constraints for subset inclusion of x+GP in {e|H.<h}
    Inputs: 
        model: Gurobi optimization model
        G: n * n_g generator matrix
        Pi:  n_pi*n matrix where {x|Pi x< 1} is the primitive polytope
        {x| Hx<h} is the set constraint
        x is the point
    Output:
        no direct output. Adds constraints to the model. 
        FUTURE: we may add lambda positive matrix as an output for debugging
    """
    (n,n_g)=G.shape
    (n_p,n)=Pi.shape
    (n_h,n)=H.shape
    Lambda=np.empty((n_h,n_p),dtype='object')
    for row in range(n_h):
        for column in range(n_p):
            Lambda[row,column]=model.addVar(lb=0)
    model.update()
    # Lambda * Pi = H * G
    for row in range(n_h):
        for column in range(n_g):
            s_left=LinExpr()
            s_right=LinExpr()
            for k in range(n_p):
                s_left.add(Lambda[row,k]*Pi[k,column])
            for k in range(n):
                s_right.add(H[row,k]*G[k,column])
            model.addConstr(s_left==s_right)
    # Lambda * 1 <= H*
    for row in range(n_h):
        s_left=LinExpr()
        s_right=LinExpr()
        for k in range(n_p):
            s_left.add(Lambda[row,k])
        for k in range(n):
            s_right.add(H[row,k]*x[k,0])
        model.addConstr(s_left<=h[row,0]-s_right)

def valuation(x):
    """
    Description: given a set of Gurobi variables, output a similar object with values
    Input:
        x: numpy vector, each val a Gurobi variable
        output: x_n: numpy array with the values 
    """
    x_n=np.empty((x.shape))
    (n_r,n_c)=x.shape
    for row in range(n_r):
        for column in range(n_c):
            x_n[row,column]=x[row,column].X
    return x_n

def pi_cube(T):
    """
        Input: T, a positive integer
        Output: Matrix Pi, characterzing unit cube
    """
    P=np.zeros((2*T,T))
    for k in range(T):
        P[2*k,k]=1
        P[2*k+1,k]=-1
    return P

def vertices_cube(T):
    from itertools import product 
    v=list(product(*zip([-1]*T,[1]*T)))
    return np.array(v)

def convexhull_RCI(s):
    """
    input: system
    dependency: execute RCI function before
    output: vertices that represent the convex hull of its computed RCI set
    """
    from scipy.spatial import ConvexHull
    v=s.x_bar.T+np.dot(s.phi,vertices_cube(s.phi.shape[1]).T).T
    return v[ConvexHull(v).vertices,:]

def visualize_RCI(s,xmin=-1,xmax=1,ymin=-1,ymax=1):
    from matplotlib.collections import PatchCollection
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import random
    if s.n!=2:
        raise(ValueError("Visualization is only for 2D state systems. Your system's state has dim=%d"%s.n))
    ax1 = plt.subplot(111)
    plt.axis('scaled')
    plt.figure(figsize=(10,10))
    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    v=convexhull_RCI(s)
    p=[patches.Polygon(v, True)]
    p=PatchCollection(p, alpha=0.4,color=(random.random(),random.random(),random.random()))
    ax1.add_collection(p)
    ax1.grid(color=(0,0,0), linestyle='--', linewidth=0.3)
    plt.show()