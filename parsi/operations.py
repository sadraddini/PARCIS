"""
@author: kasra
"""
import warnings
import numpy as np
try:
    import pypolycontain as pp
except:
    warnings.warn("You don't have pypolycontain not properly installed.")
try:
    import pydrake.solvers.mathematicalprogram as MP
    import pydrake.solvers.gurobi as Gurobi_drake
    # use Gurobi solver
    global gurobi_solver, license
    gurobi_solver=Gurobi_drake.GurobiSolver()
    license = gurobi_solver.AcquireLicense()
except:
    warnings.warn("You don't have pydrake installed properly. Methods that rely on optimization may fail.")
try:
    import parsi
except:
    warnings.warn("You don't have parsi not properly installed.")

def rci(system,order_max=10,size='min',obj='include_center'):
    """
    Given a LTI system, this function returns a rci set and its action set.
    """
    n= system.A.shape[0]
    
    for order in np.arange(1, order_max, 1/n):
        print('order',order)
        prog=MP.MathematicalProgram()
        var=parsi.rci_constraints(prog,system,T_order=order)
        T=var['T'].flatten()

        #Defining the objective function
        #objective function for minimizing the size of the RCI set
        if size=='min':
            prog.AddQuadraticErrorCost(
            Q=np.eye(len(T)),
            x_desired=np.zeros(len(T)),
            vars=T)
        #objective function for minimizing the distance between the RCI set and the set of admissible states
        if obj=='include_center':
            prog.AddQuadraticErrorCost(
            Q=np.eye(len(var['T_x'])),
            x_desired=system.X.x,
            vars=var['T_x'])

        #Result
        result=gurobi_solver.Solve(prog)    
        print('result',result.is_success())
        beta = system.beta if not 'beta' in var else result.GetSolution(var['beta'])
        if result.is_success():
            T_x=result.GetSolution(var['T_x'])
            M_x=result.GetSolution(var['M_x'])
            omega=pp.zonotope(G=result.GetSolution(var['T'])/(1-beta),x=T_x)
            theta=pp.zonotope(G=result.GetSolution(var['M'])/(1-beta),x=M_x)
            return omega,theta
        else:
            del prog
            continue    
    print("Infeasible:We couldn't find any RCI set. You can change the order_max or system.beta and try again.")
    return None,None


def viable_limited_time(system,order_max=10,size='min',obj='include_center',algorithm='slow'):
    """
    Given a LTV system, this function returns a limited time vaibale set and its action set.
    """
    number_of_steps=len(system.A)
    n= system.A[0].shape[0]
    
    for order in np.arange(1, order_max, 1/n):
        print('order',order)
        prog=MP.MathematicalProgram()
        var=parsi.viable_constraints(prog,system,T_order=order,algorithm=algorithm)
        T=[var['T'][i].flatten() for i in range(number_of_steps)]

        #Defining the objective function
        #objective function for minimizing the size of the RCI set
        if size=='min':
            [prog.AddQuadraticErrorCost(
            Q=np.eye(len(T[i])),
            x_desired=np.zeros(len(T[i])),
            vars=T[i]) for i in range(number_of_steps)]
        #objective function for minimizing the distance between the RCI set and the set of admissible states
        if obj=='include_center':
            [prog.AddQuadraticErrorCost(
            Q=np.eye(len(var['T_x'][i])),
            x_desired=system.X[i].x,
            vars=var['T_x'][i]) for i in range(number_of_steps)]

        #Result
        result=gurobi_solver.Solve(prog)    
        print('result',result.is_success())
        if result.is_success():
            T_x=[result.GetSolution(var['T_x'][i]) for i in range(number_of_steps)]
            M_x=[result.GetSolution(var['M_x']) for i in range(number_of_steps-1)]
            omega=[pp.zonotope(G=result.GetSolution(var['T'][i]),x=T_x[i]) for i in range(number_of_steps)]
            theta=[pp.zonotope(G=result.GetSolution(var['M'][i]),x=M_x[i]) for i in range(number_of_steps-1)]
            return omega,theta
        else:
            del prog
            continue    
    print("Infeasible:We couldn't find any time_limited viable set. You can change the order_max and try again.")
    return None,None




#MPC: right now it just covers point convergence
def mpc(system,horizon=1,x_desired='origin'):
    """
    MPC: Model Predictive Control
    """
    landa_terminal=10                #Terminal cost coefficient
    landa_controller=1

    n = system.A.shape[0]               # Matrix A is n*n
    m = system.B.shape[1]               # Matrix B is n*m
    if x_desired=='origin':
            x_desired=np.zeros(n)

    if system.omega==None and system.theta==None:
        system.rci()
    omega,theta=system.omega , system.theta

    prog=MP.MathematicalProgram()
    x,u=parsi.mpc_constraints(prog,system,horizon=horizon,hard_constraints=False)

    #Controller is in a shape of x=T_x + T zeta, so u=M_x + M zeta
    zeta =np.array([pp.be_in_set(prog,x[:,i],omega) for i in range(horizon)]).T                #Last one does not count
    prog.AddLinearConstraint( np.equal( theta.x.reshape(-1,1)+np.dot(theta.G,zeta) , u ,dtype='object').flatten() )

    #Objective
    
    #Cost Over x
    prog.AddQuadraticErrorCost(
    Q=np.eye(n*(horizon-1)),
    x_desired=np.tile(x_desired,horizon-1),
    vars=x[:,1:-1].T.flatten())

    #Terminal Cost
    prog.AddQuadraticErrorCost(
    Q=landa_terminal*np.eye(n),
    x_desired=x_desired,
    vars=x[:,-1].flatten())

    #Energy Cost
    prog.AddQuadraticErrorCost(
    Q=landa_controller*np.eye(m*horizon),
    x_desired=np.zeros(m*horizon),
    vars=u.flatten())

    #Result
    result=gurobi_solver.Solve(prog)    

    if result.is_success():
        u_mpc=result.GetSolution(u[:,0])
        return u_mpc


    