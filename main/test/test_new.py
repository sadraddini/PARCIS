#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:19:07 2019

@author: sadra
"""
import numpy as np

from parcis.main.system import system_LTI
from parcis.main.parcis import RCI,RCI_controller

from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes as visZ
from pypolycontain.visualization.visualize_2D import visualize_2D_zonotopes_ax as visZax

from pypolycontain.lib.zonotope import zonotope


dt=0.05
A=np.array([[1,dt],[0,1]])
B=np.array([[0,dt]]).T
W=np.array([[1,0],[0,1]])*dt*0.2
s=system_LTI(A,B,W)

s.X=zonotope(np.zeros((2,1)),  np.eye(2) *2  )
s.U=zonotope(np.zeros((1,1)),  np.eye(1) *3  )

alpha=0.5
phi,theta=RCI(s,q=10,alpha=alpha)
visZ([zonotope(np.zeros((2,1)),s.phi,color=(1,0,0))],a=0.02,title=r"RCI set")
#visZ([zonotope(np.zeros((1,1)),theta*1/(1-alpha))])

# Simulate
T=200
x,u={},{}
x[0]=np.array([-0.08,0.3]).reshape(2,1)
for t in range(T):
    print t,x[t].T
    u[t]=RCI_controller(s,x[t])
    zeta_w=np.random.random((2,1))*2-1
    x[t+1]=np.dot(A,x[t])+np.dot(B,u[t])+np.dot(W,zeta_w)

"""
Visualization
"""
import matplotlib.pyplot as plt
fig,ax = plt.subplots()
fig.set_size_inches(10, 8)
ax.set_xlabel(r"$x_0$",fontsize=20)
ax.set_ylabel(r"$x_1$",fontsize=20)
visZax(ax,[zonotope(np.zeros((2,1)),s.phi,color=(1,0,0))],a=0.02,title=r"RCI set")
ax.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)],color=(0.6,0.9,0.8),linewidth=3)
ax.plot([x[t][0,0] for t in range(T+1)],[x[t][1,0] for t in range(T+1)],'o',color=(0,0,1),markersize=7)
ax.set_title(r"",fontsize=20)