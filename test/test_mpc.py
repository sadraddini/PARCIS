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

A=np.array([[1,1],[0,1]])
B=np.array([[0],[1]])
W=pp.zonotope(G=np.eye(2),x=[0,0])*0.1
X=pp.zonotope(G=np.eye(2),x=[0,0],color='red')
U=pp.zonotope(G=np.eye(1),x=[0])*4

#Defining the system
sys=parsi.Linear_system(A,B,W=W,X=X,U=U)
sys.beta=0.9
sys.E=False
sys.rci()
sys.state=parsi.sample(sys.omega)
pp.visualize([sys.X,sys.omega])
path=sys.state.reshape(-1,1)
for step in range(40):
    #Finding the controller
    u=parsi.mpc(sys,horizon=1,x_desired='origin')
    state= sys.simulate(u)
    path=np.concatenate((path,state.reshape(-1,1)) ,axis=1)
    print(state)
    plt.plot(path[0,:],path[1,:],color='b')
    plt.pause(0.09)
plt.show()  

