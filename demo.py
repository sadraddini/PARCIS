"""
Demo:

Find a RCI set for a 2D double integrator

The safe set is unit-length box centered at origin
The set of admissible additive disturbances is 0.2-length box centered at origin
The set of admissible control inputs is 0.2-length box centered at origin

1) We find the RCI set 
2) We simulate closed-loop system starting from origin with random disturbances
"""


from rci_family import *

s=system()
s.A={}
s.B={}
s.F={}
s.g={}
s.H={}
s.r={}
s.P={}
s.q={}

s.A[1,1]=1
s.A[1,2]=1
s.A[2,1]=0
s.A[2,2]=1

s.B[1,1]=0
s.B[2,1]=1

s.F[1,1]=1
s.F[1,2]=0
s.F[2,1]=-1
s.F[2,2]=0
s.F[3,1]=0
s.F[3,2]=1
s.F[4,1]=0
s.F[4,2]=-1

s.g[1]=0.2
s.g[2]=0.2
s.g[3]=0.2
s.g[4]=0.2

s.n=2 # number of variables
s.m=1 # number of controls
s.K=10 # Design variable
s.nW=4 # Number of dis set rows
s.nX=4 # rows of X, rows of H
s.nU=2 # rows of U, rows of P

s.H=s.F
s.r=s.g
s.P[1,1]=1
s.P[2,1]=-1
s.q[1]=1
s.q[2]=1

for i in range(1,s.n+1):
	s.mu[i]=0

for j in range(1,s.m+1):
	s.v[j]=0
				
s.RCI() # Computes a set

"""
Simulation with T time steps
    - randomized disturbances
    - write state values to "state.txt"
    - initial condition = mu (centroid, default=origin)
"""
T=20
f=open("state.txt","w")
s.x=s.mu
for t in range(0,T+1):
	for i in range(1,s.n+1):
		f.write("%0.2f "%s.x[i])
	f.write("\n")	
	s.RCI_control(s.x)
	s.evolve()
f.close()
	