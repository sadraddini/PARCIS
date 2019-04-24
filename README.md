# PARCIS: Parameterized Robust Control Invariant Sets

## Introduction
Given a discrete time affine system:
<img src="http://www.sciweavers.org/tex2img.php?eq=x%5E%2B%20%3D%20Ax%20%2B%20Bu%20%2B%20Ew&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="x^+ = Ax + Bu + Ew" width="160" height="18" />
where:
* x is state, wished to be restricted to a safe polyhedron \mathbb{X}=\{x \in \mathbb{R}^n| H_x x \le h\},
* u$$ is control, with admissible set being a polyhedron \mathbb{U}=\{u \in \mathbb{R}^m| H_u u \le h_u\},
* w is additive disturbance belonging to a polytope \mathbb{W}=\{w \in \mathbb{R}^n| H_w w \le h_w\}.

The tool computes a robust control invariant set in $$\mathbb{X}$$. 


## Dependencies
* [Pypolycontain package](https://github.com/sadraddini/pypolycontain)
* Gurobi version 8.0 or later.



## Reference
The tool is partially based on the method in the following paper:
[Raković SV, Kerrigan EC, Mayne DQ, Kouramas KI. Optimized robust control invariance for linear discrete-time systems: Theoretical foundations. Automatica. 2007 May 31;43(5):831-41.](https://www.sciencedirect.com/science/article/pii/S0005109807000234):

