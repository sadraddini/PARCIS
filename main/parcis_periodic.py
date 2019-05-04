#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:39:39 2019

@author: sadra
"""

import numpy as np
from gurobipy import Model,LinExpr,QuadExpr,GRB

from pypolycontain.lib.zonotope import zonotope
from pypolycontain.lib.containment_encodings import subset_generic,subset_zonotopes

def RCI_periodic(sys,q=0,alpha=0,K=5):
    """
    Computes a Robust Control Invariant (RCI) set for Linear Periodic Systems
    """
    model=Model("RCI_periodic")
    q=K*sys.n if q==0 else q
    n_w=sys.n_w
    T=sys.T+1
    n=sys.n
    m=sys.m
    list_of_q=list(q+np.array(range(T))*n_w)
    print list_of_q
    G,theta={},{}
    for t in range(T):
        _q=list_of_q[t]
        G[t]=model.addVars(range(n),range(_q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="G_%d"%t)
        theta[t]=model.addVars(range(T),range(_q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="theta_%s"%t)
    _q=list_of_q[T-1]+n_w
    G[T]=model.addVars(range(n),range(_q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="G_%d"%T)
    model.update()
    for t in range(T):
        print "adding constraints of t",t
        A,B,W=sys.A[t],sys.B[t],sys.W[t]
        _q=list_of_q[t]
        for i in range(n):
            for j in range(_q):
                expr_x=LinExpr([(A[i,k],G[t][k,j]) for k in range(n)])
                expr_u=LinExpr([(B[i,k],theta[t][k,j]) for k in range(m)])
                model.addConstr(G[t+1][i,j]==expr_x+expr_u)
            for j in range(_q,_q+n_w):
                model.addConstr(G[t+1][i,j]==W[i,j-_q])
        G_t=np.array([G[t][i,j] for i in range(n) for j in range(_q)]).reshape(n,_q)
        theta_t=np.array([theta[t][i,j] for i in range(m) for j in range(_q)]).reshape(m,_q)
        X_t=zonotope(np.zeros((n,1)),G_t)
        U_t=zonotope(np.zeros((m,1)),theta_t)
        subset_generic(model,X_t,sys.X)
        subset_generic(model,U_t,sys.U)

    _q=list_of_q[T-1]+n_w
    for i in range(n):
        for j in range(q):
            model.addConstr(G[T][i,j+n_w*(T)]==G[0][i,j])
        for j in range(0,n_w*(T)):
            model.addConstr(G[T][i,j]==0)
    # Cost function
    J=LinExpr([(1.0/(t+1.1),G[t][i,i]) for t in range(T+1) for i in range(n)])
    model.setObjective(J)
    model.write("peridoic RCI.lp")
    model.setParam('TimeLimit', 150)
    model.optimize()
    G_num,theta_num={},{}
    for t in range(T):
        _q=list_of_q[t]
        G_num[t]=np.array([[G[t][i,j].X] for i in range(n) for j in range(_q)]).reshape(n,_q)
    _q=list_of_q[t]+n_w
    G_num[T]=np.array([[G[T][i,j].X] for i in range(n) for j in range(_q)]).reshape(n,_q)
    for t in range(T):
        _q=list_of_q[t]
        theta_num[t]=np.array([[theta[t][i,j].X] for i in range(m) for j in range(_q)]).reshape(m,_q)
    return (G_num,theta_num)


    q=K*sys.n if q==0 else q
    model=Model("RCI_periodic")
    phi=tupledict_to_array(model.addVars(range(sys.n),range(q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="phi"))
    theta=tupledict_to_array(model.addVars(range(sys.m),range(q),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="theta"))
    psi=tupledict_to_array(model.addVars(range(sys.n),range(sys.w),lb=-GRB.INFINITY,ub=GRB.INFINITY,name="psi"))
    model.update()
    _fat_matrix=np.hstack((psi,phi))
    constraints_list_of_tuples(model,[(sys.A,phi),(sys.B,theta),(-np.eye(sys.n),_fat_matrix[:,0:q])],"=")
    constraints_list_of_tuples(model,[(np.eye(sys.n),sys.W),(-np.eye(sys.n),_fat_matrix[:,q:q+sys.w])],"=")
    _outer=zonotope(np.zeros((sys.n,1)),sys.W*alpha)
    _inner=zonotope(np.zeros((sys.n,1)),psi)
    subset_generic(model,_inner,_outer)
    if sys.X!=None:
        sys.X.G*=(1-alpha)
        subset_generic(model,zonotope(np.zeros((sys.n,1)),phi),sys.X)
    if sys.U!=None:
        sys.U.G*=(1-alpha)
        subset_generic(model,zonotope(np.zeros((sys.m,1)),theta),sys.U)
    model.optimize()
    phi_n=np.array([[phi[i,j].X for i in range(phi.shape[0])] for j in range(phi.shape[1])]).T
    theta_n=np.array([[theta[i,j].X for i in range(theta.shape[0])] for j in range(theta.shape[1])]).T
#    psi_n=np.array([[psi[i,j].X for i in range(psi.shape[0])] for j in range(psi.shape[1])]).T
    sys.phi=phi_n*1.0/(1-alpha)
    sys.theta=theta_n*1.0/(1-alpha)
    return phi_n,theta_n
    
def RCI_controller(sys,x):
    """
    Based on zonotopes
    """
    model=Model("Controller")
    q=sys.phi.shape[1]
    zeta=tupledict_to_array(model.addVars(range(q),[0],lb=-1,ub=1,name="zeta"))
    zeta_abs=tupledict_to_array(model.addVars(range(q),[0],lb=0,ub=1,name="zeta_abs",obj=1))
    model.update()
    constraints_list_of_tuples(model,[(-np.eye(sys.n),x),(sys.phi,zeta)])
    model.addConstrs(zeta_abs[i,0]>=zeta[i,0] for i in range(q))
    model.addConstrs(zeta_abs[i,0]>=-zeta[i,0] for i in range(q))
    model.setParam('OutputFlag', False)
    model.optimize()
    zeta_n=np.array([zeta[i,0].X for i in range(q)]).reshape(q,1)
    print zeta_n.T
    return np.dot(sys.theta,zeta_n)
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