# PARCIS
Given a discrete time linear system

    x^+ = Ax + Bu + w
    
where 

    - x is state, wished to be restricted to a safe polyhedron S

    - u is control, with admissible set being a polyhedron

    - w is additive disturbance belonging to a polytope

We wish to find a find a robust control invariant set in S and a control policy that the state is always restricted to S. 
    
This tool is based on:

    * Raković SV, Kerrigan EC, Mayne DQ, Kouramas KI. Optimized robust control invariance for linear discrete-time systems: Theoretical foundations. Automatica. 2007 May 31;43(5):831-41.
