/*
 * Heat PDE solver using Crank Nicolson method with early exercise and entry-by-entry sor.
 *
 * Yongyi Ye
 */

#include<algorithm>
#include<eigen3/Eigen/Dense>
#include"heatpdesolver.hpp"

using namespace Eigen;
using namespace std;



MatrixXd projected_sor(MatrixXd b, double alpha);
/* Projected entry-by-entry SOR iterative method only for use in the Crank-Nicolson Heat PDE solver with
 *      early exercise. The omega, tolerance and initial guess are hard coded as well as the entries of
 *      tridiagonal matrix A in Crank-Nicolson.
 *
 * input: b: matrix b such that A * u_(m+1) = b_(m+1), used to calcule u(x, tau) at time m+1, m = 0:M-1
 *      alpha: the Courant constant, on which the entries of the tridiagonal matrix A depends.
 *
 * output: values for u(x, tau) at time m+1
 */


EarlyExCrankNicolson::EarlyExCrankNicolson(double xleft_, double xright_, double taufinal_,
                                            const Gleft &gleft_, const Gright &gright_, const Ftau &f_,
                                            const CheckEarlyExercise &checker_):
                                        EarlyExerciseSolver(xleft_, xright_, taufinal_, gleft_, gright_, f_, checker_){}

EarlyExCrankNicolson::EarlyExCrankNicolson(const EarlyExCrankNicolson &input): EarlyExerciseSolver(input){}

EarlyExCrankNicolson::~EarlyExCrankNicolson(){}

EarlyExCrankNicolson& EarlyExCrankNicolson::operator= (const EarlyExCrankNicolson &input){
    EarlyExerciseSolver::operator=(input);
    return *this;
}

MatrixXd EarlyExCrankNicolson::solve_pde(int n, int m){
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

    // Set up the matrix B - A will be set up in the projected_sor function, it is
    //      a tridiagonal matrix depends only on alpha.
    MatrixXd B(MatrixXd::Zero(n-1, n-1));
    for(int row = 0; row < n-1; row++){ B(row, row) = 1 - alpha; }
    for(int row = 0; row < n-2; row++){ B(row, row+1) = 0.5 * alpha; }
    for(int row = 1; row < n-1; row++){ B(row, row-1) = 0.5 * alpha; }

    // Using Crank Nicolson method to compute the nodes.
    MatrixXd b(n-1, 1);
    MatrixXd u_next(n-1, 1);

    for(int row = 1; row < m+1; row++){

      for(int i = 1; i < n; i++){ b(i-1, 0) = u(row-1, i); }
      b =  B * b;
      b(0, 0) = b(0, 0) +  0.5 * alpha * (u(row, 0) + u(row-1,0));
      b(n-2, 0) = b(n-2, 0) + 0.5 * alpha * (u(row, n) + u(row-1, n));

      u_next = projected_sor(b, alpha);

      // Check for early exercise and copy entries to the mesh
      for(int i = 1; i < n; i++){ u(row, i) = max<double>(u_next(i-1, 0), checker->check_early_exercise(n, m, i, row)); }
    }

    return u;
}

HeatPdeSolver* EarlyExCrankNicolson::clone() const{
    return new EarlyExCrankNicolson(*this);
}


/********** projected sor ***********/
MatrixXd projected_sor(MatrixXd b, double alpha){
    //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
    //no need for the M, N matrix in this case - use entries of A directly

    //get the dimension of A - same as the number of rows of b
    int n = b.rows();


    /** the hard-coded parameters - change if necessary! **/
    // initial guess
    //MatrixXd x0(MatrixXd::Constant(n,1,1));
    MatrixXd x0(MatrixXd::Constant(n,1,0));
    // omega
    double w = 1.2;
    // tolerance
    double tol = 0.000001;
    /**    end    **/


    //set up matrix A
    MatrixXd A(MatrixXd::Zero(n-1, n-1));
    for(int row = 0; row < n-1; row++){ A(row, row) = 1 + alpha; }
    for(int row = 0; row < n-2; row++){ A(row, row+1) = -0.5 * alpha; }
    for(int row = 1; row < n-1; row++){ A(row, row-1) = -0.5 * alpha; }

    //stop criteria
    MatrixXd r0(b - A*x0);
    double stop_iter_resid = tol * r0.norm();
    MatrixXd r(r0);

    //counter
    int ic = 0;

    //the two vectors in the iteration
    MatrixXd x_old(MatrixXd::Constant(n,1, 7));
    MatrixXd x_new(x0);

    //start the iteration
    while( (x_new-x_old).norm()>tol && ic <= INT_MAX-1){    //consecutive approx. stopping
    //while( r.norm()>stop_iter_resid && ic <= INT_MAX-1){  //residual-based stopping
        x_old = x_new;

        //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
        for(int j=0; j<n; j++){
            //update he entries of x_new one by one - band width of A is 2
            double x_new_j = 0;
            for(int k = max(j-1,0); k<j; k++){ x_new_j =  x_new_j + A(j,k)*x_new(k,0); }
            for(int k = j+1; k<min(j+2,n); k++){ x_new_j =  x_new_j + A(j,k)*x_old(k,0); }
            x_new(j,0) = (1-w)*x_old(j,0) + w*( -x_new_j + b(j,0))/A(j,j);
        }

        r = b - A*x_new;
        ic ++;
    }

    //if(ic == INT_MAX-1){
        //cout << "SOR does not converge" << endl;
    //}
    //cout << endl << "sor op count" << ic << endl;
    return x_new;

}
