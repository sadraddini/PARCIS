#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:40:15 2019

@author: sadra
"""
import numpy as np
from gurobipy import Model,LinExpr,QuadExpr,GRB

from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.containment_encodings import subset_generic,subset_zonotopes

def RCI(sys,q=0,alpha=0,K=5):
    """
    Computes a Robust Control Invariant (RCI) set for LTI system sys
    """
    q=K*sys.n if q==0 else q
    model=Model("RCI")
    phi=tupledict_to_array(model.addVars(range(sys.n),range(q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="phi"))
    theta=tupledict_to_array(model.addVars(range(sys.m),range(q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="theta"))
    psi=tupledict_to_array(model.addVars(range(sys.n),range(sys.w),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="psi"))
#    alpha=tupledict_to_array(model.addVars([0],[0],lb=0,ub=1,name="alpha"))
    model.update()
    _fat_matrix=np.hstack((psi,phi))
    constraints_list_of_tuples(model,[(sys.A,phi),(sys.B,theta),(-np.eye(sys.n),_fat_matrix[:,0:q])],"=")
    constraints_list_of_tuples(model,[(np.eye(sys.n),sys.W),(-np.eye(sys.n),_fat_matrix[:,q:q+sys.w])],"=")
    _outer=zonotope(np.zeros((sys.n,1)),sys.W*alpha)
    _inner=zonotope(np.zeros((sys.n,1)),psi)
#    subset_generic(model,_inner,_outer)
    subset_zonotopes(model,_inner,_outer)
    model.optimize()
#    print phi
    phi_n=np.array([[phi[i,j].X for i in range(phi.shape[0])] for j in range(phi.shape[1])]).T
    theta_n=np.array([[theta[i,j].X for i in range(theta.shape[0])] for j in range(theta.shape[1])]).T
    psi_n=np.array([[psi[i,j].X for i in range(psi.shape[0])] for j in range(psi.shape[1])]).T
    _fat_matrix_n=np.array([[_fat_matrix[i,j].X for i in range(_fat_matrix.shape[0])] for j in range(_fat_matrix.shape[1])]).T
    print "fat_matrix \n",_fat_matrix_n
    print "hstack of psi and phi \n",np.hstack((psi_n,phi_n))
    return phi_n,theta_n,psi_n
    
    
"""
Auxilary Gurobi Shortcut Functions
used for convenience
"""

def tupledict_to_array(mytupledict):
    # It should be 2D
    n,m=max(mytupledict.keys())
    n+=1
    m+=1
    array=np.empty((n,m),dtype="object")
    for i in range(n):
        for j in range(m):
            array[i,j]=mytupledict[i,j]
    return array

def constraints_list_of_tuples(model,mylist,sign="="):
    term_0=mylist[0]
    ROWS,COLUMNS=term_0[0].shape[0],term_0[1].shape[1]
    for row in range(ROWS):
        for column in range(COLUMNS):
            expr=LinExpr()
            for term in mylist:
                q,qp=term[0].shape[1],term[1].shape[0]
                if q!=qp:
                    raise ValueError(term,"q=%d qp=%d"%(q,qp))
                if type(term[1][0,column])==type(model.addVar()):
                    expr.add(LinExpr([(term[0][row,k],term[1][k,column]) for k in range(q)]))
                elif type(term[0][row,0])==type(model.addVar()):
                    expr.add(LinExpr([(term[1][k,column],term[0][row,k]) for k in range(q)]))
                else:
                    expr.addConstant(sum([term[1][k,column]*term[0][row,k] for k in range(q)]))
            if sign=="<":
                model.addConstr(expr<=0)
            elif sign=="=":
                model.addConstr(expr==0)
            elif sign==">=":
                model.addConstr(expr>=0)
            else:
                raise "sign indefinite"