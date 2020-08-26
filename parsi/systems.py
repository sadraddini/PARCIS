"""
@author: kasra
Defining classes for linear systems both LTV (linear time variant) and LTI (linear time invariant)
"""
import numpy as np
try:
    import pypolycontain as pp
except:
    print('Error: pypolycontain package is not installed correctly') 


class Linear_system:
    
    def __init__(self,A,B,W=None,X=None,U=None):
        """
        A and B are the matrices for the dynamics of the linear syst: x^+ = Ax+Bu+w
        W is the disturbance set: w \in W
        X is the state space set: x \in X
        U is the control input set: u \in U
        This class covers both LTI and LTV systems: if teh system is LTV the inputs need to be in the form of lists, otherwise they need to be numpy.
        """
        if type(A)==list and type(B)==list:
            self.sys = 'LTV'                #it's a linear time-variant system

        elif type(A)!=list and type(B)!=list:
            if W == None:
                self.sys = 'non-disturbed'              # For x^+ = Ax+Bu (without given disturbance)
            else:
                self.sys = 'LTI'                #it's a linear time-invariant system
                assert(len(A) == W.x.shape[0]), "The dimension of matrix A should be the same as dimension of W"
                assert(len(B) == W.x.shape[0]), "The dimension of matrix B should be the same as dimension of W"

            assert(len(A) == len(B)), "The number of rows in matrix A and B need to be the same"

        else:
            raise ValueError ('Input areguments need to be in the form of a list for LTV systems and np.ndarray for LTI systems')

        self.A = np.array(A)
        self.B = np.array(B)
        self.W = W
        self.X = X
        self.U = U
        self.state = np.array([0,0])
        self.beta=0.2 if self.sys=='LTI' else None               #For finding RCI
        self.E=True if self.sys=='LTI' else None             #For finding RCI
        self.omega=None             #RCI set
        self.theta=None             #Action set

    def __repr__(self):
        return self.sys

    def simulate(self,u):
        x_next = np.dot(self.A,self.state)+np.dot(self.B,u)+sample(self.W)
        self.state = x_next
        return x_next

    def rci(self,order_max=10,size='min',obj='include_center'):
        import parsi
        omega,theta=parsi.rci(self,order_max=10,size=size,obj=obj)
        self.omega=omega
        self.theta=theta
        return omega,theta


#sampling from a set represented in zonotope
def sample(zonotope):
    random_box=np.random.uniform(low=-1.0, high=1.0, size=zonotope.G.shape[1])
    output= zonotope.x + np.dot(zonotope.G,random_box)
    return output


# class Coupled_linear_sytem(Linear_system):
#     """
#     This class covers coupled linear systems in the form:
#     x^+_i = A_{ii}x_i + B_{ii}u_i + w_i + A_{ij}x_j + B_{ij}u_j
#     """
#     def __init__(self,system,partition):




        