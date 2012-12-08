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

    // Set up the boundary conditions for the (m+1)*(n+1) mesh
    MatrixXd u(m+1, n+1);

    for(int i = 0; i < n+1; i++){ u(0, i) = (*f)(xleft + i*delta_x); }
    for(int j = 1; j < m+1; j++){
      u(j, 0) = (*gleft)(j * delta_t);
      u(j, n) = (*gright)(j * delta_t);
    }

    // Set up the matrix A.
    MatrixXd A(MatrixXd::Zero(n-1, n-1));

    for(int row = 0; row < n-1; row++){ A(row, row) = 1 + 2*alpha; }
    for(int row = 0; row < n-2; row++){ A(row, row+1) = -alpha; }
    for(int row = 1; row < n-1; row++){ A(row, row-1) = -alpha; }

    // Using backward Euler to compute the nodes.
    MatrixXd b(n-1, 1);
    MatrixXd u_next(n-1, 1);

    for(int row = 1; row < m+1; row++){

      for(int i = 1; i < n; i++){ b(i-1, 0) = u(row-1, i); }

      b(0, 0) = b(0, 0) + alpha * u(row, 0);
      b(n-2, 0) = b(n-2, 0) + alpha * u(row, n);

      u_next = solver->solve(A, b);

      for(int j = 1; j < n; j++){ u(row, j) = u_next(j-1, 0); }
    }

    return u;
}

HeatPdeSolver* BackwardEuler::clone() const{
    return new BackwardEuler(*this);
}
