#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:19:07 2019

@author: sadra
"""
import numpy as np

from parcis.main.system import system_LTI
from parcis.main.parcis import RCI

dt=0.03
A=np.array([[1,dt],[0,1]])
B=np.array([[0,dt]]).T
W=np.eye(2)*0.1
s=system_LTI(A,B,W)
phi,theta,psi=RCI(s,q=4,alpha=0.5)