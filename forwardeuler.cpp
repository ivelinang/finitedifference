/*
 * Heat PDE solver using Forward Euler method
 *
 * Yongyi Ye
 */

#include"heatpdesolver.hpp"

#include<stdlib.h>
using namespace std;

#include<eigen3/Eigen/Dense>
using namespace Eigen;

ForwardEuler::ForwardEuler(double xleft_, double xright_, double taufinal_,
                               const Gleft &gleft_, const Gright &gright_, const Ftau &f_):
                            HeatPdeSolver(xleft_, xright_, taufinal_, gleft_, gright_, f_){}

ForwardEuler::ForwardEuler(const ForwardEuler &input):HeatPdeSolver(input){}

ForwardEuler::~ForwardEuler(){}

ForwardEuler& ForwardEuler::operator= (const ForwardEuler &input){
    HeatPdeSolver::operator=(input);
}

MatrixXd ForwardEuler::solve_pde(int n, int m){
    double delta_t = taufinal / static_cast<double>(m);
    double delta_x = (xright - xleft) / static_cast<double>(n);
    double alpha = delta_t / (delta_x * delta_x);

    // Set up the boundary conditions.
    MatrixXd u;
    VectorXd left_bc, right_bc;

    u.setZero(m + 1, n - 1);
    left_bc.setZero(m + 1);
    right_bc.setZero(m + 1);

    for(int i = 0; i < n - 1; ++i){ u(0, i) = (*f)(xleft + (i+1) * delta_x); }
    for(int j = 0; j < m + 1; ++j){
      left_bc(j) = (*gleft)(j * delta_t);
      right_bc(j) = (*gright)(j * delta_t);
    }

    // Using Forward Euler to compute the nodes.
    VectorXd b, u_next;

    for(int row = 0; row < m; row++){
      u(row+1, 0) = alpha * u(row, 1) + (1 - 2*alpha) * u(row, 0) + alpha * left_bc(row);

      for(int col = 1; col < n - 2; col++){
          u(row+1, col) = alpha * u(row, col+1) + (1 - 2*alpha) * u(row, col) + alpha * u(row, col-1);
      }

      u(row+1, n-2) = alpha * right_bc(row) + (1 - 2*alpha) * u(row, n-2) + alpha * u(row, n-3);
    }

    return u;
}

