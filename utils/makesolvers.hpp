/*
 * Declaration.
 *
 * This script is a collection of make functions that produce and return the Black-Scholes
 * partial differential equation solvers (../blackscholes/bspdesolve.hpp) for either
 * American or European put and call given the relevant parameters (stock price, strike,
 * volatility, maturity, interest, continuous dividend rate).
 *
 * Yongyi Ye
 */

#ifndef MAKESOLVERS_HPP
#define MAKESOLVERS_HPP

#include"../blackscholes/bspdesolver.hpp"

/* Note:
 *      euro - European
 *      amer - American
 *      call - Call Option
 *       put - Put Option
 *
 *        fe - Forward Euler
 *        be - Backward Euler
 *        cn - Crank Nicolson
 *
 *  cholesky - solve the linear system created by be/cn at each time step using Cholesky decomposition & solve
 *       sor - solve the linear system created by be/cn at each time step using Successive Over Relaxation
 *                  iterative method with given omega (w, default 1.2) and tolerance (tol, default 10^-6).
 */


/* European Call */
BsPdeSolver make_euro_call_fe(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_euro_call_be_lu(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_euro_call_be_cholesky(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_euro_call_be_sor(double s, double k, double vol, double t, double r, double q,
                                    double w = 1.2, double tol = 0.000001);

BsPdeSolver make_euro_call_cn_cholesky(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_euro_call_cn_sor(double s, double k, double vol, double t, double r, double q,
                                    double w = 1.2, double tol = 0.000001);


/* European Put */
BsPdeSolver make_euro_put_fe(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_euro_put_be_cholesky(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_euro_put_be_sor(double s, double k, double vol, double t, double r, double q,
                                    double w = 1.2, double tol = 0.000001);

BsPdeSolver make_euro_put_cn_cholesky(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_euro_put_cn_sor(double s, double k, double vol, double t, double r, double q,
                                    double w = 1.2, double tol = 0.000001);


/* American Call */
BsPdeSolver make_amer_call_fe(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_amer_call_cn(double s, double k, double vol, double t, double r, double q,
                                double w = 1.2, double tol = 0.000001);


/* American Put */
BsPdeSolver make_amer_put_fe(double s, double k, double vol, double t, double r, double q);

BsPdeSolver make_amer_put_cn(double s, double k, double vol, double t, double r, double q,
                                double w = 1.2, double tol = 0.000001);

#endif
