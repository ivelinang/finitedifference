/*
 * Implementation of the Black-Scholes PDE solver for European options
 *
 * Yongyi Ye
 */


#include<math.h>
#include<iostream>
#include<eigen3/Eigen/Dense>
#include"bspdesolver.hpp"

using namespace std;
using namespace Eigen;

BsPdeSolver::BsPdeSolver(double s_, double k_, double vol_, double t_, double r_, double q_, const HeatPdeSolver &solver_):
                            s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_), solver(solver_){

    //compute value for a and b
    a = (r-q) / (vol*vol) - 0.5;
    double b_ = (r-q) / (vol*vol) + 0.5;
    b = pow(b_, 2) + 2*q / (vol*vol);
}

MatrixXd BsPdeSolver::solve_pde(double alpha, int M){
    int N = compute_n(alpha, M);
    MatrixXd soln(solver->solve_pde(N, M)); // soln is a (M+1)*(N+1) matrix
    return soln;
}

double BsPdeSolver::compute_price(double alpha, int M){
    MatrixXd soln(this->solve_pde(alpha, M));
    int N = compute_n(alpha, M);

    // find i such that x_i <= x_compute <= x_i+1
    double x_compute = log(s/k);
    double x_left = log(s/k) + (r - q - vol*vol/2.0) * t - 3 * vol * sqrt(t);
    // x_right - x_left = 6 * vol * sqrt(t) = N * delta_x
    double delta_x = (6 * vol * sqrt(t)) / static_cast<double>(N);
    // truncate the floating part of computed i
    int i = (x_compute - x_left)/delta_x;

    double tau_final = t * vol * vol / 2.0;

    if(i+1 > N){ return -1; }

    double x_1 = x_left + i * delta_x;
    double x_2 = x_left + (i+1) * delta_x;

    double v_1 = exp(-a*x_1 - b*tau_final) * soln(M, i);
    double v_2 = exp(-a*x_2 - b*tau_final) * soln(M, i+1);

    double s_1 = k * exp(x_1);
    double s_2 = k * exp(x_2);

    /*
    cout << "x_compute " << x_compute << endl
         << "x_i " << x_1 << endl << "x_ii " << x_2 << endl
         << "s_i " << s_1 << endl << "s_ii " << s_2 << endl
         << "v_i " << v_1 << endl << "v_ii " << v_2 << endl;
    */

    double price = ((s_2 - s) * v_1 + (s - s_1) * v_2) / (s_2 - s_1);
    return price;
}

/* forward finite difference approx. of delta */
double BsPdeSolver::compute_delta(double alpha, int M){
    MatrixXd soln(this->solve_pde(alpha, M));
    int N = compute_n(alpha, M);

    // find i such that x_i <= x_compute <= x_i+1
    double x_compute = log(s/k);
    double x_left = log(s/k) + (r - q - vol*vol/2.0) * t - 3 * vol * sqrt(t);
    // x_right - x_left = 6 * vol * sqrt(t) = N * delta_x
    double delta_x = (6 * vol * sqrt(t)) / static_cast<double>(N);
    // truncate the floating part of computed i
    int i = (x_compute - x_left)/delta_x;

    double tau_final = t * vol * vol / 2.0;

    if(i+1 > N){ return -1; }

    double x_1 = x_left + i * delta_x;
    double x_2 = x_left + (i+1) * delta_x;

    double v_1 = exp(-a*x_1 - b*tau_final) * soln(M, i);
    double v_2 = exp(-a*x_2 - b*tau_final) * soln(M, i+1);

    double s_1 = k * exp(x_1);
    double s_2 = k * exp(x_2);

    double delta = (v_2 - v_1) / (s_2 - s_1);
    return delta;
}

/* forward finite difference approx. of gamma */
double BsPdeSolver::compute_gamma(double alpha, int M){
    MatrixXd soln(this->solve_pde(alpha, M));
    int N = compute_n(alpha, M);

    // find i such that x_i <= x_compute <= x_i+1
    double x_compute = log(s/k);
    double x_left = log(s/k) + (r - q - vol*vol/2.0) * t - 3 * vol * sqrt(t);
    // x_right - x_left = 6 * vol * sqrt(t) = N * delta_x
    double delta_x = (6 * vol * sqrt(t)) / static_cast<double>(N);
    // truncate the floating part of computed i
    int i = (x_compute - x_left)/delta_x;

    double tau_final = t * vol * vol / 2.0;

    if(i+2 > N || i-1 < 0){ return -1; }

    double x_0 = x_left + (i-1) * delta_x;
    double x_1 = x_left + i * delta_x;
    double x_2 = x_left + (i+1) * delta_x;
    double x_3 = x_left + (i+2) * delta_x;

    double v_0 = exp(-a*x_0 - b*tau_final) * soln(M, i-1);
    double v_1 = exp(-a*x_1 - b*tau_final) * soln(M, i);
    double v_2 = exp(-a*x_2 - b*tau_final) * soln(M, i+1);
    double v_3 = exp(-a*x_3 - b*tau_final) * soln(M, i+2);

    double s_0 = k * exp(x_0);
    double s_1 = k * exp(x_1);
    double s_2 = k * exp(x_2);
    double s_3 = k * exp(x_3);

    double gamma = ( (v_3 - v_2) / (s_3 - s_2) - (v_1 - v_0) / (s_1 - s_0) ) /
                   ( (s_3 + s_2) / 2.0 - (s_1 + s_0) / 2.0);

    return gamma;
}


/* forward finite difference approx. of theta */
double BsPdeSolver::compute_theta(double alpha, int M){
    MatrixXd soln(this->solve_pde(alpha, M));
    int N = soln.cols();

    // find i such that x_i <= x_compute <= x_i+1
    double x_compute = log(s/k);
    double x_left = log(s/k) + (r - q - vol*vol/2.0) * t - 3 * vol * sqrt(t);
    // x_right - x_left = 6 * vol * sqrt(t) = N * delta_x
    double delta_x = (6 * vol * sqrt(t)) / static_cast<double>(N);
    // truncate the floating part of computed i
    int i = (x_compute - x_left)/delta_x;

    double tau_final = t * vol * vol / 2.0;
    double delta_tau = tau_final/static_cast<double>(M);

    if(i+1 > N){ return -1; }

    double x_1 = x_left + i * delta_x;
    double x_2 = x_left + (i+1) * delta_x;

    double v_1 = exp(-a*x_1 - b*tau_final) * soln(M, i);
    double v_2 = exp(-a*x_2 - b*tau_final) * soln(M, i+1);
    double v_3 = exp(-a*x_1 - b * (tau_final-delta_tau)) * soln(M-1, i);
    double v_4 = exp(-a*x_2 - b * (tau_final-delta_tau)) * soln(M-1, i+1);

    double s_1 = k * exp(x_1);
    double s_2 = k * exp(x_2);

    double price_1 = ((s_2 - s) * v_1 + (s - s_1) * v_2) / (s_2 - s_1);
    double price_2 = ((s_2 - s) * v_3 + (s - s_1) * v_4) / (s_2 - s_1);

    double theta = (price_2 - price_1)/delta_tau;
    return theta;
}

int BsPdeSolver::compute_n(double alpha, int M){
    double tau_final = t * vol * vol / 2.0;
    double delta_tau = tau_final/static_cast<double>(M);

    // x_right - x_left = 6 * vol * sqrt(t)
    double rl = 6 * vol * sqrt(t);

    double float_n = rl / sqrt(delta_tau/alpha);
    int n = static_cast<int>(float_n);
    return n;
}




