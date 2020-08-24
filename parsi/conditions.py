"""
@author: kasra
"""
import numpy as np
try:
    import pypolycontain as pp
except:
    print('Error: pypolycontain package is not installed correctly') 


def rci_constraints(program,system,T_order=3):
    """
    This function adds necessary constraints for RCI (robust control invariant) set for a LTI.
    """
    assert(system.sys=='LTI'), "The system has to be LTI (linear time invariant)"
    n = system.A.shape[0]               # Matrix A is n*n
    m = system.B.shape[1]               # Matrix B is n*m
    k = round(n*T_order)
    p = system.W.G.shape[1]
    
    E=system.E
    flag=0
    if E==True:
        e=program.NewContinuousVariables( n, p, 'E')                #E is n*p
        if (system.X.x==0).all() and (system.U.x==0).all() and (system.W.x==0).all():
            flag=1
    elif E==False:
        e=np.zeros((n,p))
        system.beta=0
    else:
        raise ValueError('E in the argument needs to be either False or True')      
    beta = system.beta

    #Defining Variables
    T=program.NewContinuousVariables( n, k, 'T')                #T is n*k
    M=program.NewContinuousVariables( m, k, 'M')                #M is m*k
    T_x=program.NewContinuousVariables( n,'T_center')               #T_x is the x_bar in the paper
    M_x=program.NewContinuousVariables( m,'M_center')               #M_x is the u_bar in the paper

    #Defining Constraints
    left_side = np.concatenate( (np.dot(system.A,T)+np.dot(system.B,M) , system.W.G) ,axis=1)               #[AT+BM,W]
    right_side = np.concatenate( (e,T),axis=1)              #[E,T]
    program.AddLinearConstraint(np.equal(left_side,right_side,dtype='object').flatten())             #[AT+BM,W]==[E,T]

    if flag==1:
        _,_,beta = pp.subset(program, pp.zonotope(G=e,x=np.zeros(e.shape[0])) , pp.zonotope(G=system.W.G,x=np.zeros(system.W.G.shape[0])) , alpha='scalar')             #Z(0,E) \subset Z(0,beta*W)
        program.AddBoundingBoxConstraint(0,1,beta)              #CHECK IT, I NEED 0<=beta<1
        # program.AddLinearConstraint(beta < 1)
        # program.AddLinearConstraint(beta >= 0)
        program.AddLinearConstraint(np.equal(T_x, np.zeros(T.shape[0]),dtype='object').flatten())
        program.AddLinearConstraint(np.equal(M_x, np.zeros(M.shape[0]),dtype='object').flatten())
        if system.X!=None:
            _,_,beta1=pp.subset(program, pp.zonotope(G=T,x=np.zeros(T.shape[0])) , system.X , alpha='scalar' )
            program.AddLinearConstraint(np.equal(beta1,(1-beta),dtype='object'))
        if system.U!=None:
            _,_,beta2=pp.subset(program, pp.zonotope(G=M,x=np.zeros(M.shape[0])) , system.U , alpha='scalar' )
            program.AddLinearConstraint(np.equal(beta2,(1-beta),dtype='object'))
    elif flag==0:
        if E==True:
            pp.subset(program, pp.zonotope(G=e,x=np.zeros(e.shape[0])) , pp.zonotope(G=beta *system.W.G,x=np.zeros(system.W.G.shape[0])) )
        if system.X!=None:
            pp.subset(program, pp.zonotope(G=T/(1-beta),x=T_x) , system.X )
        if system.U!=None:
            pp.subset(program, pp.zonotope(G=M/(1-beta),x=M_x) , system.U )
        center=np.equal( np.dot(system.A , T_x) + np.dot(system.B , M_x) + system.W.x , T_x ,dtype='object').flatten()
        program.AddLinearConstraint(center)               #A*x_bar+ B*u_bar +w_bar==x_bar       #IT WILL MAKE PROBLEM WHEN IT BECOMES TRUE!

    #print('cccccccccccenter',center)
    #print('ceeeeeeeeeenter',center==True)
    #center=np.delete(center , center==True)
    
    output={
        'T':T,
        'T_x': T_x,
        'M':M,
        'M_x': M_x
    }
    if flag==1:
        output['beta']=beta
    
    return output


def mpc_constraints(program,system,horizon=1,hard_constraints=True):

    if system.sys=='LTI':
        n = system.A.shape[0]               # Matrix A is n*n
        m = system.B.shape[1]               # Matrix B is n*m

        x=program.NewContinuousVariables( n,horizon,'x')
        x=np.concatenate((system.state.reshape(-1,1),x ), axis=1)
        u=program.NewContinuousVariables( m,horizon,'u')

        #Dynamics
        dynamics=np.equal( np.dot(system.A , x[:,:-1]) + np.dot(system.B , u) , x[:,1:] ,dtype='object').flatten()
        program.AddLinearConstraint(dynamics)

        #Hard constraints over state
        if hard_constraints==True:
            if system.X!=None:
                _=[pp.be_in_set(program,x[:,i],system.X) for i in range(1,horizon)]
            #Hard constraints over control input
            if system.U!=None:
                _=[pp.be_in_set(program,u[:,i],system.U) for i in range(horizon)]

        return x,u