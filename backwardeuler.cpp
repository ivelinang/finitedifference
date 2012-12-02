/*
 * Heat PDE solver using Backward Euler method
 *
 * Yongyi Ye
 */

#include"heatpdesolver.hpp"

#include<stdlib.h>
using namespace std;

#include<eigen3/Eigen/Dense>
using namespace Eigen;

BackwardEuler::BackwardEuler(double xleft_, double xright_, double taufinal_,
                               const Gleft &gleft_, const Gright &gright_, const Ftau &f_,
                               const LinearSolver &solver_):
                            HeatPdeSolver(xleft_, xright_, taufinal_, gleft_, gright_, f_),
                               solver(solver_){}

BackwardEuler::BackwardEuler(const BackwardEuler &input):HeatPdeSolver(input),
                                                         solver(input.solver){}

BackwardEuler::~BackwardEuler(){}

BackwardEuler& BackwardEuler::operator= (const BackwardEuler &input){
    HeatPdeSolver::operator=(input);
    solver = input.solver;
}

MatrixXd BackwardEuler::solve_pde(int n, int m){
    double delta_t = taufinal / static_cast<double>(m);
    double delta_x = (xright - xleft) / static_cast<double>(n);
    double alpha = delta_t / (delta_x * delta_x);

    // Set up the boundary conditions.
    MatrixXd u;
    VectorXd left_bc, right_bc;

    u.setZero(m+1, n-1);
    left_bc.setZero(m+1);
    right_bc.setZero(m+1);

    for(int i = 0; i < n-1; ++i){ u(0, i) = (*f)(xleft + (i+1) * delta_x); }
    for(int j = 0; j < m+1; ++j){
      left_bc(j) = (*gleft)(j * delta_t);
      right_bc(j) = (*gright)(j * delta_t);
    }

    // Set up the matrix A.
    MatrixXd A;
    A.setZero(n-1, n-1);

    for(int row = 0; row < n-1; row++){ A(row, row) = 1 + 2*alpha; }
    for(int row = 0; row < n-2; row++){ A(row, row+1) = -alpha; }
    for(int row = 1; row < n-1; row++){ A(row, row-1) = -alpha; }

    // Using forward Euler to compute the nodes.
    VectorXd b, u_next;

    for(int row = 1; row < m+1; row++){
      b.setZero(n-1);
      u_next.setZero(n-1);

      for(int i = 0; i < n-1; i++){ b(i) = u(row-1, i); }

      b(0) = b(0) + alpha * left_bc(row);
      b(n-2) = b(n-2) + alpha * right_bc(row);

      u_next = solver->solve(A, b);

      for(int j = 0; j < n - 1; ++j)
        u(row, j) = u_next(j);
    }

    return u;
}

