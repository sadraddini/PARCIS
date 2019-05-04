#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 18:41:43 2019

@author: sadra
"""

class system_LTI:
    def __init__(self,A,B,W):
        """
        Linear discrete-time system
        """
        self.A,self.B,self.W=A,B,W
        self.n,self.m=self.B.shape
        self.w=W.shape[1]
        if A.shape!=(self.n,self.n) or W.shape[0]!=self.n:
            raise ValueError 
        self.X=None
        self.U=None
            
    def set_X(self,X):
        self.X=X
    
    def set_U(self,U):
        self.U
        
class system_T_periodic:
    def __init__(self,A,B,W):
        """
        Linear discrete-time system
        """
        self.A,self.B,self.W=A,B,W
        self.n,self.m=self.B[0].shape
        self.n_w=W[0].shape[1]
        self.X=None
        self.U=None
        self.T=max(self.A.keys() )
            
    def set_X(self,X):
        self.X=X
    
    def set_U(self,U):
        self.U
        
