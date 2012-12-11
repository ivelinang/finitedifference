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
    return *this;
}

MatrixXd ForwardEuler::solve_pde(int n, int m){
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

    // Using Forward Euler to compute the nodes.
    for(int row = 1; row < m+1; row++){
      for(int col = 1; col < n; col++){
          u(row, col) = alpha * u(row-1, col+1) + (1 - 2*alpha) * u(row-1, col) + alpha * u(row-1, col-1);
      }
    }

    return u;
}


HeatPdeSolver* ForwardEuler::clone() const{
    return new ForwardEuler(*this);
}
