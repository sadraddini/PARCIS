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
            Q=np.eye(len(var['T_x'])) if size=='min' else -1*np.eye(len(T)),
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

    return "Infeasible:We couldn't find any RCI set. You can change the order_max or system.beta and try again."