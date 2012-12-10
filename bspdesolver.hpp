/*
 * Declaration.
 *
 * The solver to solve the Black Scholes partial differential equation to get the price for
 * an European/American call or put given the parameters of the option and a heat pde solver
 * of choice (forward euler, backward euler, or crank-nicolson).
 *
 * The heat pde is in the form of u_tau = u_xx such that u is a function of x and tau with
 * boundary conditions:
 * -- f(x) for tau = 0, x_left < x < x_right
 * -- gleft(tau) for x = x_left, 0 < tau < tau_final
 * -- gright(tau) for x = x_right, 0< tau < tau_final
 *
 * This solver contains a function that takes in the courant constant (alpha) and number of
 * partitions M on the tau axis to compute the number of partitions on the x axis N, in
 * order to form a mesh of M*N to use finite difference method of solve the heat pde.
 *
 * Yongyi Ye
 */

#ifndef BLACKSCHOLESPDE_HPP
#define BLACKSCHOLESPDE_HPP

#include"wrapper.hpp"
#include"heatpdesolver.hpp"

class BsPdeSolver{
/* the solver is for European options AND American Options depends on the heat pde solver passed in */

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

        Wrapper<HeatPdeSolver> solver;  // the heat pde solver

        // constants involved to transform a black scholes pde into a heat pde
        double a;
        double b;

    public:
        /* constructor */
        BsPdeSolver(double s_, double k_, double vol_, double t_, double r_, double q_, const HeatPdeSolver &solver_);

        /* Function to solve the pde and compute the option value, delta, gamma, or theta at time 0 given the
         *      asset price at time 0 - s and the courant constant alpha_temp, mesh size M on the tau axis.
         * It will start by figuring out the mesh size N on the x axis and call the pde solver. */
        MatrixXd solve_pde(double alpha, int M);

        double compute_price(double alpha, int M);
        double compute_delta(double alpha, int M);
        double compute_gamma(double alpha, int M);
        double compute_theta(double alpha, int M);

        /* Function to be used by the solve_pde to figure out N */
        int compute_n(double alpha, int M);

};

#endif
