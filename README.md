# PARCIS: Parameterized Robust Control Invariant Sets

## Introduction
Given a discrete time affine system:
$$x^+ = Ax + Bu + c + w$$
where:
* $$x$$ is state, wished to be restricted to a safe polyhedron $$\mathbb{X}=\{x \in \mathbb{R}^n| H_x x \le h\}$$,
* $$u$$ is control, with admissible set being a polyhedron $$\mathbb{U}=\{u \in \mathbb{R}^m| H_u u \le h_u\}$$,
* $$w$$ is additive disturbance belonging to a polytope $$\mathbb{W}=\{w \in \mathbb{R}^n| H_w w \le h_w\}$$.

The tool computes a robust control invariant set in $$\mathbb{X}$$. 


## Dependencies
Gurobi version 6.0 or later. Python 3.5 or later.


## What's new in version 1.1
1. Matrices are numpy-based, which results in slightly faster pre-processing. 

## Reference
The tool is partially based on the method in the following paper:
[Raković SV, Kerrigan EC, Mayne DQ, Kouramas KI. Optimized robust control invariance for linear discrete-time systems: Theoretical foundations. Automatica. 2007 May 31;43(5):831-41.](https://www.sciencedirect.com/science/article/pii/S0005109807000234):

