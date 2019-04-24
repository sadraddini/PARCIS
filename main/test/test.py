#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 15:35:19 2018

@author: sadra
"""
from ana_PARCIS import system,convexhull_RCI,np,visualize_RCI

s=system(2,1)

s.A=np.array([[1,1],[0,1]])
s.B=np.array([[0,1]]).T
s.c=np.array([[0,0]]).T

s.H=np.array([[1,0],[0,1],[-1,0],[0,-1]])
s.h=np.array([[1,1,1,1]]).T

s.F=np.array([[1,-1]]).T
s.f=np.array([[1,1]]).T

s.W=np.array([[2,0],[0,3]])*(1/10)
s.RCI()

visualize_RCI(s)