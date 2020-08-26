"""
@author: kasra
"""
import numpy as np
import matplotlib.pyplot as plt
try:
    import pypolycontain as pp
except:
    print('Error: pypolycontain package is not installed correctly') 
try:
    import parsi
except:
    raise ModuleNotFoundError("parsi package is not installed properly")

number_of_steps=15
A=[np.array([[1+0.001*t,1],[0,1-0.001*t]]) for t in range(number_of_steps)]
B=[np.array([[0],[1+0.002*t]]) for t in range(number_of_steps)]
W=[pp.zonotope(G=np.eye(2),x=[0,0])*0.2 for t in range(number_of_steps)]
X=[pp.zonotope(G=np.eye(2)*3,x=[0,0],color='red') for t in range(number_of_steps)]
U=[pp.zonotope(G=np.eye(1),x=[0]) for t in range(number_of_steps)]

sys=parsi.Linear_system(A,B,W=W,X=X,U=U)

omega,theta=parsi.viable_limited_time(sys,order_max=10,algorithm='slow')

shift=10
for i in range(number_of_steps):
    omega[i].x=omega[i].x+10*i*np.array([1,0])
    X[i].x=X[i].x+10*i*np.array([1,0])
pp.visualize([*X,*omega])
plt.show()