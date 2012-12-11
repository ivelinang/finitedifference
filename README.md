finitedifference
================

A simple finite difference implementation to solve the Black-Scholes PDE

The solver to solve the Black Scholes partial differential equation to get the price for
an European/American call or put given the parameters of the option and a heat pde solver
of choice (forward euler, backward euler, or crank-nicolson for European, 
special forward euler, special crank-nicolson with entry-by-entry sor for American).

The heat pde is in the form of u_tau = u_xx such that u is a function of x and tau with
boundary conditions:
-- f(x) for tau = 0, x_left < x < x_right
-- gleft(tau) for x = x_left, 0 < tau < tau_final
-- gright(tau) for x = x_right, 0< tau < tau_final

This solver contains a function that takes in the courant constant (alpha) and number of
partitions M on the tau axis to compute the number of partitions on the x axis N, in
order to form a mesh of MxN to use finite difference method of solve the heat pde.

__Yongyi Ye__


*Structure*

blackscholes: the folder contains the closed-form Black-Scholes formula calculator (bs.hpp/cpp); the solver to solve
              Black-Scholes PDE with the market parameters and a solver of choice( bspdesolver.hpp/cpp); a checker for checking
              early exercise when computing values for American options (checkearlyexercise.hpp/cpp); and the gleft, gright, f
              classes that defines boundary conditions for the heat pde.

finitedifference: the folder of the main function and makefile.

heatpdesolvers: this folder contains a file for all Heat PDE solver's declarations (heatpdesolver.hpp), 
                and the implementation for each specific solver.(mentioned above)

linearalgebra: this folder contains the solver and underlying utilities to be used to solve for x in a linear system of
               Ax = b with various methods. (all declarations are in linearsolver.hpp)

utils: contains a class that deals with pointers and memory management (new/delete) using generic programming.


*Timeline:*

Nov.30: Transformed previously written linear solver methods into classes with a common interface
        that has a function called solve to solve for x in linear system Ax = b.

Dec.1: Made sure the Linear Solver classes are working properly and correctly.
       Finished building the heat PDE solver structure. (Forward Euler, Backward Euler, Crank-Nicolson)

Dec.2: Finished the implementation for the core heat pde solver. Made sure the numbers are correct.
       Fixed index errors and made code more robust by creating less objects.

Dec.7: Added Black-Scholes PDE solver using the Heat PDE (u_tau = u_xx) solver and the boundary conditions
       for European put and call options.

Dec.9: Added boundary conditions, early exercise checker, special Heat PDE solver for solving Heat PDE with
       early exercise features for American put and call options. The Black-Scholes PDE solver can be used
       to solve for American options by passing in an American-style-specific Heat PDE solver.


       Things to be done next: performance improvement and make function for different solvers.
