/*
 * Heat PDE solver using Crank-Nicolson method
 *
 * Yongyi Ye
 */

#include"heatpdesolver.hpp"

#include<stdlib.h>
using namespace std;

#include<eigen3/Eigen/Dense>
using namespace Eigen;

CrankNicolson::CrankNicolson(double xleft_, double xright_, double taufinal_,
                               const Gleft &gleft_, const Gright &gright_, const Ftau &f_,
                               const LinearSolver &solver_):
                            HeatPdeSolver(xleft_, xright_, taufinal_, gleft_, gright_, f_),
                               solver(solver_){}

CrankNicolson::CrankNicolson(const CrankNicolson &input):HeatPdeSolver(input),
                                                         solver(input.solver){}

CrankNicolson::~CrankNicolson(){}

CrankNicolson& CrankNicolson::operator= (const CrankNicolson &input){
    HeatPdeSolver::operator=(input);
    solver = input.solver;
}

MatrixXd CrankNicolson::solve_pde(int n, int m){
    double delta_t = taufinal / static_cast<double>(m);
    double delta_x = (xright - xleft) / static_cast<double>(n);
    double alpha = delta_t / (delta_x * delta_x);

    // Set up the boundary conditions.
    MatrixXd u;
    VectorXd left_bc, right_bc;

    u.setZero(m + 1, n - 1);
    left_bc.setZero(m + 1);
    right_bc.setZero(m + 1);

    for(int i = 0; i < n - 1; ++i)
      u(0, i) = (*f)(xleft + (i + 1) * delta_x);
    for(int j = 0; j < m + 1; ++j){
      left_bc(j) = (*gleft)(j * delta_t);
      right_bc(j) = (*gright)(j * delta_t);
    }

    // Set up the matrix A and B.
    MatrixXd A, B;
    A.setZero(n - 1, m - 1);
    B.setZero(n - 1, m - 1);

    for(int row = 0; row < n; ++row)
      A(row, row) = 1 + alpha;
    for(int row = 0; row < n - 1; ++row)
      A(row, row + 1) = - 0.5 * alpha;
    for(int row = 1; row < n; ++row)
      A(row, row - 1) = - 0.5 * alpha;

    // Using Crank Nicolson method to compute the nodes.
    VectorXd b, u_next;

    for(int row = 0; row < m; ++row){
      b.setZero(n - 1);
      u_next.setZero(n - 1);

      for(int i = 0; i < n - 1; ++i)
        b(i) = u(row - 1, i);

      b *= B;

      b(0) += 0.5 * alpha * left_bc(row) + 0.5 * alpha * left_bc(row + 1);
      b(n - 2) += 0.5 * alpha * right_bc(row) + 0.5 * alpha * right_bc(row + 1);

      u_next = solver->solve(A, b);

      for(int j = 0; j < n - 1; ++j)
        u(row + 1, j) = u_next(j);
    }

    return u;
}

