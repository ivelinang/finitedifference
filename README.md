finitedifference
================

A simple finite difference implementation to solve the Black-Scholes PDE

*Timeline:*
Nov.30: Transformed previously written linear solver methods into classes with a common interface
        that has a function called solve to solve for x in linear system Ax = b.

Dec.1: Made sure the Linear Solver classes are working properly and correctly.
       Finished building the heat PDE solver structure. (Forward Euler, Backward Euler, Crank-Nicolson)
