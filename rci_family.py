""" 
Version 1.0

Synthesizing Robust Control Invariant (RCI) Sets for linear 
discrete-time systems subject to additive disturbances

Based on paper:
"Rakovic, S. V., Eric C. Kerrigan, David Q. Mayne, and Konstantinos I. Kouramas.
"Optimized robust control invariance for linear discrete-time systems:
Theoretical foundations." Automatica 43, no. 5 (2007): 831-841. "

Developed by:
	Sadra Sadraddini
	Boston University, Boston, MA
	sadra at bu dot edu
	February 2017


Input:	Given matrix A, B and polytope W (Fw<g),
		state set X (Hx<r), control set U (Pv<q),
		
Goal:	find a RCI set \Omega \subseteq X, and a 
		invariance inducing control stratgey \Omega -> U
		
Output: Return matrices in M and D, parameters of RCI set and the controller 

Notes:
....... Usage with proper citation is allowed
....... The code is still under development
"""

from gurobipy import *
import random

class system:
	def __init__(self):
		self.A={} #	 A matrix, enter as a dictionary
		self.B={} #	 B matrix, enter as a dictionary
		self.F={} #	 F, Fw<=g
		self.g={} #	 
		
		self.H={} #	 H x-H mu<beta*r-delta
		self.mu={} #  
		self.r={}
		self.beta=1
		
		self.P={} #	 Pu - Pv < gamma*q-eps
		self.v={} #	 
		self.q={}
		self.gamma=1
		
		self.alpha=0.0 # Contraction Factor
		
		self.n=0 # number of variables
		self.m=0 # number of controls
		self.K=1 # Design variable, degree
		self.nW=1 # Number of dis set rows
		self.nX=1 # rows of X, rows of H
		self.nU=1 # rows of U, rows of P
		self.AA={} # AA[k,i,j] is i,j entry of A^k
		self.HAB={} # HAB[k,i,j] is i,j entry of H*A^k*B
		self.FAB={} # FAB[k,i,j] is i,j entry of F*A^k*B
		self.M={} 
		self.D={}
		self.xbar={}
		self.ubar={}
		self.x={}
		self.u={}
		
	
	def compute_HAB(self):
		for k in range(0,self.K):
			for i in range(1,self.nX+1):
				for j in range(1,self.m+1):
					self.HAB[k,i,j]=0
					for q in range(1,self.n+1):
						for p in range(1,self.n+1):
							self.HAB[k,i,j]+=self.H[i,q]*self.AA[k,q,p]*self.B[p,j]
	def compute_FAB(self):
		for k in range(0,self.K):
			for i in range(1,self.nW+1):
				for j in range(1,self.m+1):
					self.FAB[k,i,j]=0
					for q in range(1,self.n+1):
						for p in range(1,self.n+1):
							self.FAB[k,i,j]+=self.F[i,q]*self.AA[k,q,p]*self.B[p,j]
						
	
	def compute_AA(self):
		for i in range(1,self.n+1):
			for j in range(1,self.n+1):
				self.AA[0,i,j]=int(i==j)
				self.AA[1,i,j]=self.A[i,j]
		for k in range(0,self.K):
			for i in range(1,self.n+1):
				for j in range(1,self.n+1):
					self.AA[k+1,i,j]=0
					for p in range(1,self.n+1):
						self.AA[k+1,i,j]+=self.A[i,p]*self.AA[k,p,j]
						
				
								
	def compute_D(self):
		for k in range(0,self.K):
			for i in range(1,self.n+1):
				for j in range(1,self.n+1):
					self.D[k,i,j]=self.AA[self.K-k-1,i,j]
					for p in range(1,self.n+1):
						for q in range(1,self.m+1):
							for eta in range(0,self.K-k-1):
								self.D[k,i,j]+=self.AA[self.K-k-2-eta,i,p]*self.B[p,q]*self.M[eta,q,j]
							
					  
		
		
	def RCI(self):
		M={}
		GAM={}
		PI={}
		LAMBDA={}
		xbar={}
		ubar={}
		model=Model("RCIS")
		
		# These are new!
		self.beta=model.addVar(lb=0,ub=GRB.INFINITY)
		self.gamma=model.addVar(lb=0,ub=GRB.INFINITY)
		self.rho=model.addVar(lb=0,ub=GRB.INFINITY,obj=1)
		
		
		for i in range(1,self.n+1):
			xbar[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
		for i in range(1,self.m+1):
			ubar[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
			
		for k in range(0,self.K):
			for i in range(1,self.m+1):
				for j in range(1,self.n+1):
					M[k,i,j]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
		for i in range(1,self.nW+1):
			for j in range(1,self.nW+1):
				LAMBDA[i,j]=model.addVar(lb=0)
				
		for k in range(0,self.K):
			for i in range(1,self.nX+1):
				for j in range(1,self.nW+1):
					GAM[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)
			for i in range(1,self.nU+1):
				for j in range(1,self.nW+1):
					PI[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)		
			
		model.update()
		
		# A xbar + B u bar= xbar
		for i in range(1,self.n+1):
			s=LinExpr()
			for j in range(1,self.n+1):
				s.addTerms(self.A[i,j], xbar[j])
			for j in range(1,self.m+1):
				s.addTerms(self.B[i,j], ubar[j])
			model.addConstr(s == xbar[i])
		
		for i in range(1,self.n+1):
			model.addConstr(xbar[i]==0)
		for j in range(1,self.m+1):
			model.addConstr(ubar[j]==0)
		
		# Equation 4.9-1
		# Lambda * g <= alpha * g
		for i in range(1,self.nW+1):
			s=LinExpr()
			for j in range(1,self.nW+1):
				s.addTerms(self.g[j],LAMBDA[i,j])
			model.addConstr(s <= self.alpha*self.g[i])
		
		# Equation 4.9-2
		for i in range(1,self.nW+1):
			for j in range(1,self.n+1):
				leftHand=LinExpr()
				for p in range(1,self.nW+1):
					leftHand.addTerms(self.F[p,j],LAMBDA[i,p])
				# Now the heavy part
				rightHand=LinExpr()
				for eta in range(0,self.K):
					for p in range(1,self.m+1): 
						rightHand.addTerms(self.FAB[self.K-1-eta,i,p],M[eta,p,j])
				righthandConstant=0
				for p in range(1,self.n+1):
					righthandConstant+=self.F[i,p]*self.AA[self.K,p,j]
				model.addConstr( leftHand == rightHand + righthandConstant)
		
		# Equation 4.10a)
		for i in range(1,self.nX+1):
			Hmu=0
			for j in range(1,self.n+1):
				Hmu+=self.H[i,j]*self.mu[j]
			s=LinExpr()
			Hx=LinExpr()
			for j in range(1,self.n+1):
				Hx.addTerms(-self.H[i,j],xbar[j])
			for k in range(0,self.K):
				for j in range(1,self.nW+1):
					s.addTerms(self.g[j],GAM[k,i,j])
			model.addConstr(s <= (1-self.alpha)*self.beta*self.r[i] + (1-self.alpha)*Hmu + (1-self.alpha)*Hx )
		  
		# Equation 4.10b)
		for k in range(0,self.K):
			for i in range(1,self.nX+1):
				for j in range(1,self.n+1):
					leftHand=LinExpr()
					for p in range(1,self.nW+1):
						leftHand.addTerms(self.F[p,j],GAM[k,i,p])
					# Now the heavy part
					rightHand=LinExpr()
					for eta in range(0,k):
						for p in range(1,self.m+1): 
							rightHand.addTerms(self.HAB[k-1-eta,i,p],M[eta,p,j])
					righthandConstant=0
					for p in range(1,self.n+1):
						righthandConstant+=self.H[i,p]*self.AA[k,p,j]
					model.addConstr( leftHand == rightHand + righthandConstant)
		
		# Equation 4.11a)
		for i in range(1,self.nU+1):
			Pv=0
			for j in range(1,self.m+1):
				Pv+=self.P[i,j]*self.v[j]
			s=LinExpr()
			Pu=LinExpr()
			for j in range(1,self.m+1):
				Pu.addTerms(-self.P[i,j],ubar[j])
			for k in range(0,self.K):
				for j in range(1,self.nW+1):
					s.addTerms(self.g[j],PI[k,i,j])
			model.addConstr(s <= (1-self.alpha)*self.gamma*self.q[i]+ (1-self.alpha)*Pv + (1-self.alpha)*Pu )
		
		# Equation 4.11b)
		for k in range(0,self.K):
			for i in range(1,self.nU+1):
				for j in range(1,self.n+1):
					leftHand=LinExpr()
					for p in range(1,self.nW+1):
						leftHand.addTerms(self.F[p,j],PI[k,i,p])
					# Now the heavy part
					rightHand=LinExpr()
					for p in range(1,self.m+1): 
						rightHand.addTerms(self.P[i,p],M[k,p,j])
					model.addConstr( leftHand == rightHand)			   
		
		model.addConstr( self.rho >= self.beta)
		model.addConstr( self.rho >= self.gamma)	
		
		model.optimize()
		if model.Status==3:
			return False
		else:
			print "xbar is:"
			for i in range(1,self.n+1):
				print "-",i,":",xbar[i].X,"-",
			print ""
			for i in range(1,self.n+1):
				self.xbar[i]=xbar[i].X
			for j in range(1,self.m+1):
				self.ubar[j]=ubar[j].X
			for k in range(0,self.K):
				for i in range(1,self.m+1):
					for j in range(1,self.n+1):
						self.M[k,i,j]=M[k,i,j].X
			self.beta=self.beta.X
			self.gamma=self.gamma.X
			return True
					
	
	
	
	
	
	
					
	######################### CHECK ##########################	 
	######################### CHECK ##########################					
	def is_RCI(self,x):
		w={}
		model=Model("RCIS_check")
		for k in range(0,self.K):
			for i in range(1,self.n+1):
				w[k,i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
		
			
		model.update()
		# x bar = x + D * w
		for i in range(1,self.n+1):
			s=LinExpr()
			for k in range(0,self.K):
				for p in range(1,self.n+1):
					s.addTerms(self.D[k,i,p], w[k,p])
			model.addConstr(s + self.xbar[i] == x[i])
	
		
		# w in (1-alpha)^-1 W
		for k in range(0,self.K):
			for i in range(1,self.nW+1):
				s=LinExpr()
				for p in range(1,self.n+1):
					s.addTerms(self.F[i,p],w[k,p])
				model.addConstr( (1-self.alpha)*s <= self.g[i])
		
		model.optimize()
		if model.Status==3:
			return False
		else:
			return True
	######################### END ##########################				  
	######################### END ##########################				  
				 
				
				
	######################### BEGIN CONTROL ##########################	 
	def RCI_control(self,x):
		w={}
		model=Model("RCIS_check")
		for k in range(0,self.K):
			for i in range(1,self.n+1):
				w[k,i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name="w(%d,%d)" %(k,i))
				   
		model.update()
		# x = xbar + D * w
		for i in range(1,self.n+1):
			s=LinExpr()
			for k in range(0,self.K):
				for p in range(1,self.n+1):
					s.addTerms(self.D[k,i,p], w[k,p])
			model.addConstr(s + self.xbar[i] == x[i])
	
		
		# w in (1-alpha)^-1 W
		for k in range(0,self.K):
			for i in range(1,self.nW+1):
				s=LinExpr()
				for p in range(1,self.n+1):
					s.addTerms(self.F[i,p],w[k,p])
				model.addConstr( (1-self.alpha)*s <= self.g[i])
		
		J = QuadExpr()
		for k in range(0,self.K):
			for i in range(1,self.n+1):
				J.add(w[k,i]*w[k,i])
#		  print J
		model.setObjective(J,GRB.MINIMIZE)
		model.optimize()
		u={}
		for i in range(1,self.m+1):
			u[i]=self.ubar[i]
			for k in range(0,self.K):
				for p in range(1,self.n+1):
					u[i]+=self.M[self.K-1-k,i,p]*w[k,p].X
		self.u=u
		return u
	######################### END CONTROL ##########################					 
	
	######################### BEGIN IMPLEMENATION ##########################				  
	def evolve(self):
		w={}
		model=Model("RCIS_check")
		for i in range(1,self.n+1):
			w[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,obj=random.random()-0.5)
		model.update()
		for i in range(1,self.nW+1):
			s=LinExpr()
			for p in range(1,self.n+1):
				s.addTerms(self.F[i,p], w[p])
			model.addConstr( s <= self.g[i])
		model.optimize()
		print "randomized the disturbance!:",
		xnew={}
		for i in range(1,self.n+1):
			xnew[i]=w[i].X*0.99#*random.random()
			print xnew[i],
			for p in range(1,self.n+1):
				xnew[i]+=self.A[i,p]*self.x[p]
			for p in range(1,self.m+1):
				xnew[i]+=self.B[i,p]*self.u[p]
		self.x=xnew
	######################### END IMPLEMENATION ##########################					

	def RCI_vertex(self,direction):
		w={}
		x={}
		model=Model("RCIS_check")
		for k in range(0,self.K):
			for i in range(1,self.n+1):
				w[k,i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
		
		for i in range(1,self.n+1):
			x[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY, obj=direction[i])
			
		model.update()
		# x bar = x + D * w
		for i in range(1,self.n+1):
			s=LinExpr()
			for k in range(0,self.K):
				for p in range(1,self.n+1):
					s.addTerms(self.D[k,i,p], w[k,p])
			model.addConstr(s + self.xbar[i] == x[i])
		
#		model.addConstr(x[4]==0)
#		model.addConstr(x[2]==0)
		
		# w in (1-alpha)^-1 W
		for k in range(0,self.K):
			for i in range(1,self.nW+1):
				s=LinExpr()
				for p in range(1,self.n+1):
					s.addTerms(self.F[i,p],w[k,p])
				model.addConstr( (1-self.alpha)*s <= self.g[i])
		
		model.optimize()
		vertex={}
		for i in range(1,self.n+1):
			vertex[i]=x[i].X
		return vertex			
					
		
		