/*
 * Heat PDE solver using Forward Euler method with early exercise.
 *
 * Yongyi Ye
 */

#include<algorithm>
#include<eigen3/Eigen/Dense>
#include"heatpdesolver.hpp"

using namespace std;
using namespace Eigen;


EarlyExForwardEuler::EarlyExForwardEuler(double xleft_, double xright_, double taufinal_,
                                             const Gleft &gleft_, const Gright &gright_, const Ftau &f_,
                                             const CheckEarlyExercise &checker_):
                                EarlyExerciseSolver(xleft_, xright_, taufinal_, gleft_, gright_, f_, checker_){}

EarlyExForwardEuler::EarlyExForwardEuler(const EarlyExForwardEuler &input): EarlyExerciseSolver(input){}

EarlyExForwardEuler::~EarlyExForwardEuler(){}

EarlyExForwardEuler& EarlyExForwardEuler::operator= (const EarlyExForwardEuler &input){
    EarlyExerciseSolver::operator=(input);
    return *this;
}


MatrixXd EarlyExForwardEuler::solve_pde(int n, int m){
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
          double num_1 = alpha * u(row-1, col+1) + (1 - 2*alpha) * u(row-1, col) + alpha * u(row-1, col-1);
          double num_2 = checker->check_early_exercise(n, m, col, row);
          u(row, col) = max<double>(num_1, num_2);
      }
    }

    return u;
}

HeatPdeSolver* EarlyExForwardEuler::clone() const{
    return new EarlyExForwardEuler(*this);
}
