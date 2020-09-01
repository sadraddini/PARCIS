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

number_of_subsystems=2
n=2*number_of_subsystems
m=1*number_of_subsystems

A=np.random.rand(n,n)* 0.3
#B=np.random.rand(n,m)* 0.0001
#A=np.zeros((n,n))
B=np.zeros((n,m))
A[0:2,0:2]= np.array([[1,1],[0,1]]) 
A[2:4,2:4]= np.array([[1,1],[0,1]]) 
B[0:2,0]= np.array([0,1])
B[2:4,1]= np.array([0,1])

W=pp.zonotope(G=np.eye(n)*0.1,x=np.zeros(n))
X=pp.zonotope(G=np.eye(n),x=np.zeros(n),color='red')
U=pp.zonotope(G=np.eye(m),x=np.zeros(m))

system=parsi.Linear_system(A,B,W=W,X=X,U=U)
sub_sys=parsi.sub_systems(system,partition_A=[2]*number_of_subsystems,partition_B=[1]*number_of_subsystems)

for i in range(number_of_subsystems):
    sub_sys[i].U.G=np.array([sub_sys[i].U.G])

omega,theta = parsi.decentralized_rci(sub_sys)

#pp.visualize(omega)
sub_sys[0].X.color='red'
pp.visualize([sub_sys[0].X,omega[0]])
plt.show()
