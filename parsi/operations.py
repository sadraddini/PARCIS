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
            M_x=[result.GetSolution(var['M_x'][i]) for i in range(number_of_steps-1)]
            omega=[pp.zonotope(G=result.GetSolution(var['T'][i]),x=T_x[i]) for i in range(number_of_steps)]
            theta=[pp.zonotope(G=result.GetSolution(var['M'][i]),x=M_x[i]) for i in range(number_of_steps-1)]
            return omega,theta
        else:
            del prog
            continue    
    print("Infeasible:We couldn't find any time_limited viable set. You can change the order_max and try again.")
    return None,None



def sub_systems(system,partition_A,partition_B):
    """
    The input is a large system and partition over system.A 
    example: [2,3,1] creates A_{11}=2*2, A_{22}=3*3, and A_{33}=1*1
    """
    assert(len(partition_A)==len(partition_B)), "length of vector partition_A has to be equal to the length of the vector partition_B"
    from itertools import accumulate
    number_of_subsys=len(partition_A)
    par_accum_A=list(accumulate(partition_A))
    par_accum_B=list(accumulate(partition_B))
    par_accum_A.insert(0,0)
    par_accum_B.insert(0,0)

    A=[[ system.A[par_accum_A[i]:par_accum_A[i+1] , par_accum_A[j]:par_accum_A[j+1]] for i in range(number_of_subsys)] for j in range(number_of_subsys)]
    B=[[ system.B[par_accum_A[i]:par_accum_A[i+1] , par_accum_B[j]:par_accum_B[j+1]] for i in range(number_of_subsys)] for j in range(number_of_subsys)]

    X=pp.decompose(system.X,partition_A)
    U=pp.decompose(system.U,partition_B)

    #Right now, it just covers disturbances with order 1.
    W = pp.boxing_order_reduction(system.W,desired_order=1)
    W=pp.decompose(W,partition_A)

    sys=[ parsi.Linear_system(A[i][i] , B[i][i] , W=W[i],X=X[i] , U=U[i]) for i in range(number_of_subsys)]
    for i in range(number_of_subsys):
        for j in range(number_of_subsys):
            if i==j:
                continue
            if (A[i][j]==0).all()==False:
                sys[i].A_ij[j]= A[i][j]
            if (B[i][j]==0).all()==False:
                sys[i].B_ij[j]= B[i][j]

    return sys


def decentralized_rci(sys,method='centralized',initial_guess='nominal'):
    """
    Given a set of coupled linear sub-systems
    """
    for i in sys:
        i.E=False               #Setting all E=False
    order_max=10
    number_of_subsys= len(sys)

    if method=='centralized':
        
        if initial_guess=='nominal':
            X_i,U_i=[],[]
            for i in range(number_of_subsys):
                omega,theta = rci(sys[i])
                X_i.append(omega)
                U_i.append(theta)
            
            number_of_columns=[sys[i].W.G.shape[1] for i in range(number_of_subsys)]
            for i in range(number_of_subsys):
                for j in range(number_of_subsys):
                    if j in sys[i].A_ij:
                        number_of_columns[i]=number_of_columns[i]+ X_i[j].G.shape[1]
                    if j in sys[i].B_ij:
                        number_of_columns[i]=number_of_columns[i]+ U_i[j].G.shape[1]
            
        n=[len(sys[i].A) for i in range(number_of_subsys)]
        disturb=[]
        for i in range(number_of_subsys):
            disturb.append(sys[i].W)
        W_i=disturb
        print('W_iiiiiiiiiiii',W_i[0].G)

        for order in np.arange(6, order_max, 1/max(n)):
            
            prog= MP.MathematicalProgram()
            
            #Disturbance over each sub-system
            d_aug=[ prog.NewContinuousVariables( n[i], number_of_columns[i], 'd') for i in range(number_of_subsys) ]
            d_aug_x=[prog.NewContinuousVariables( n[i],'d_center') for i in range(number_of_subsys) ]

            #setting the disturbance for each sub-system
            for i in range(number_of_subsys):                
                sys[i].W = pp.zonotope(G=d_aug[i],x=d_aug_x[i])

            #rci_constraints
            var= [parsi.rci_constraints(prog,sys[i],T_order=order) for i in range(number_of_subsys)]

            #Correctness Criteria
            alpha_x = [pp.subset(prog,pp.zonotope(G= var[i]['T'], x=var[i]['T_x']), X_i[i] , alpha='vector')[2] for i in range(number_of_subsys)]
            [prog.AddBoundingBoxConstraint(0,np.inf,alpha_x[i]) for i in range(number_of_subsys)]
            alpha_u = [pp.subset(prog,pp.zonotope(G= var[i]['M'], x=var[i]['M_x']), U_i[i] , alpha='vector')[2] for i in range(number_of_subsys)]
            [prog.AddBoundingBoxConstraint(0,np.inf,alpha_u[i]) for i in range(number_of_subsys)]

            #Computing the disturbance set
            
            for i in range(number_of_subsys):
                for j in range(number_of_subsys):
                    if j in sys[i].A_ij:
                        w=pp.zonotope(G= np.dot(np.dot( sys[i].A_ij[j], X_i[j].G), np.diag(alpha_x[j]) ) , x= np.dot( sys[i].A_ij[j], X_i[j].x))
                        disturb[i] = disturb[i]+ w
                    if j in sys[i].B_ij:
                        w=pp.zonotope(G= np.dot(np.dot( sys[i].B_ij[j], U_i[j].G), np.diag(alpha_u[j]) ) , x= np.dot( sys[i].B_ij[j], U_i[j].x))
                        disturb[i] = disturb[i]+ w
            print("THE ORDER IS=",order)
            
            #Disturbance
            [prog.AddLinearConstraint(np.equal(d_aug_x[i], disturb[i].x,dtype='object').flatten()) for i in range(number_of_subsys)]
            [prog.AddLinearConstraint(np.equal(d_aug[i], disturb[i].G,dtype='object').flatten()) for i in range(number_of_subsys)]

            #Result
            result=gurobi_solver.Solve(prog)    
            print('result',result.is_success())
            if result.is_success():
                T_x=[result.GetSolution(var[i]['T_x']) for i in range(number_of_subsys)]
                M_x=[result.GetSolution(var[i]['M_x']) for i in range(number_of_subsys)]
                omega=[pp.zonotope(G=result.GetSolution(var[i]['T']),x=T_x[i]) for i in range(number_of_subsys)]
                theta=[pp.zonotope(G=result.GetSolution(var[i]['M']),x=M_x[i]) for i in range(number_of_subsys)]
                return omega,theta
            else:
                del prog
                print('W_i',W_i[0].G)
                disturb=W_i
                print("I AM HEEEEEERRRRR")
                continue    
        print("Infeasible:We couldn't find a set of decentralized rci sets. You can change the order_max and try again.")
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


    