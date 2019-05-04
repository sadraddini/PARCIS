#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 11:00:54 2019

@author: sadra
"""

import numpy as np

from parcis.main.system import system_T_periodic
from parcis.main.parcis_periodic import RCI_periodic

from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZax

from pypolycontain.lib.zonotope import zonotope

dt=0.05
A,B,W={},{},{}
A[0]=np.array([[1,dt],[0,1]])
B[0]=np.array([[0,dt]]).T
W[0]=np.array([[1,0],[0,1]])*dt*0.4

A[1]=np.array([[1,dt],[1,0.6]])
B[1]=np.array([[0,dt]]).T
W[1]=np.array([[1,0],[0,1]])*dt*0.2

A[2]=np.array([[1,dt],[-1,1]])
B[2]=np.array([[0,dt]]).T
W[2]=np.array([[1,0],[0,1]])*dt*0.4

A[3]=np.array([[1,dt],[-0.5,1.2]])
B[3]=np.array([[0,dt]]).T
W[3]=np.array([[1,0],[0,1]])*dt*0.2

s=system_T_periodic(A,B,W)

s.X=zonotope(np.zeros((2,1)),  np.eye(2) *2  )
s.U=zonotope(np.zeros((1,1)),  np.eye(1) *3  )

G,theta=RCI_periodic(s,q=12)

T=s.T+1
for t in range(T+1):
    if t==T:
        visZ([zonotope(np.zeros((2,1)),G[t],color=(t/(s.T+1.5),0,1-t/(s.T+1.5)))],a=0.02,title=r"$\mathcal{Z}(0,G_{0})=\mathcal{Z}(0,G_{%d})$"%T)
    else:
        visZ([zonotope(np.zeros((2,1)),G[t],color=(t/(s.T+1.5),0,1-t/(s.T+1.5)))],a=0.02,title=r"$\mathcal{Z}(0,G_{%d})$"%t)

