/*
 * Implementation of the linear solver classes that uses iterative method to get
 * solutions, including:
 *
 *      Jacobi
 *      JacobiBanded
 *      GaussSeidel
 *      GaussSeidelBanded
 *      SOR
 *      SOR-Banded
 *
 * they all inherit from the base class interface linearsolver.
 *
 */

#include<math.h>
#include<iostream>
#include<eigen3/Eigen/Dense>

#include"linearsolver.hpp"
#include"decomposer.hpp"

using namespace std;
using namespace Eigen;

/** Jacobi **/
Jacobi::Jacobi(double tol_): tol(tol_){}

MatrixXd Jacobi::solve(MatrixXd A, MatrixXd b){
    /*
     * function to solve x for linear equation Ax = b using Jacobi iteration method
     * converge when A is strictly diagonally dominant,  weakly diagonally dominant and irreducible, or consistently ordered
     *
     * input: A: a n*n matrix
     *        b: a n*1 vector
     *        tol: a tolerance factor to stop the iteration
     *
     * output: a n*1 vector as the solution to Ax=b
     *
     */

    //for the matrix of M and N, where M = D and N = L_A + U_A (M+N = A)
    int n = A.rows();

    MatrixXd M(MatrixXd::Zero(n,n));
    for(int i=0; i<n; i++){
        M(i,i) = A(i,i);
    }
    MatrixXd N(A-M);

    //initial guess
    MatrixXd x0(MatrixXd::Zero(n,1));
    //MatrixXd x0(MatrixXd::Constant(n,1,1));

    //stop criteria
    MatrixXd r0(b - A*x0);
    double stop_iter_resid = tol * r0.norm();


    MatrixXd x(n,1);
    x = x0;

    MatrixXd r(n,1);
    r = r0;

    MatrixXd D_inv(M.inverse());

    MatrixXd b_new(n,1);
    b_new = D_inv * b;

    //counter
    int ic = 0;

    //the two vectors in the iteration
    MatrixXd x_old(MatrixXd::Zero(n,1));
    MatrixXd x_new(x0);

    //start the iteration
    //while( (x_new-x_old).norm()>tol && ic <= INT_MAX-1){    //consecutive approx. stopping
    while( r.norm()>stop_iter_resid && ic <= INT_MAX-1){  //residual-based stopping
        x_old = x_new;
        x_new = - D_inv * N * x_old + b_new;
        //x = - D_inv * N * x + b_new;
        r = b - A*x_new;
        ic ++;

        //print out the first n iterations
        //if(ic <=3){ cout << "Jacobi " << ic << " iteration: " << endl << x_new << endl << endl; }
    }

    if(ic == INT_MAX-1){
        cout << "Jacobi does not converge" << endl;
    }

    //cout << "jacobi op count:" << ic << endl;
    return x_new;
}

LinearSolver* Jacobi::clone() const{
    return new Jacobi(*this);
}


/* Gauss Seidel */
GaussSeidel::GaussSeidel(double tol_):tol(tol_){}

MatrixXd GaussSeidel::solve(MatrixXd A, MatrixXd b){
    /*
     * function to solve x for linear equation Ax = b using gauss-siedel iteration method
     * converge when A is spd, strictly diagonally dominant, diagonally dominant and irreducible, or consistently ordered
     *
     * input: A: a n*n matrix
     *        b: a n*1 vector
     *        tol: a tolerance factor to stop the iteration
     *
     * output: a n*1 vector as the solution to Ax=b
     *
     */

    //construct the matrix D, L_A, U_A, where L_A + D + U_A = A
    int n = A.rows();

    MatrixXd D(MatrixXd::Zero(n,n));
    MatrixXd L_A(MatrixXd::Zero(n,n));
    MatrixXd U_A(MatrixXd::Zero(n,n));

    for(int i=0; i<n; i++){
        D(i,i) = A(i,i);

        for(int j=0; j<i; j++){
            L_A(i,j) = A(i,j);
        }

        for(int j=i+1; j<n; j++){
            U_A(i,j) = A(i,j);
        }
    }
    //initial guess
    MatrixXd x0(MatrixXd::Zero(n,1));
    //MatrixXd x0(MatrixXd::Constant(n,1,1));

    //stop criteria
    MatrixXd r0(b - A*x0);
    double stop_iter_resid = tol * r0.norm();

    //initialize variables for the iteration
    MatrixXd x(x0);
    MatrixXd r(r0);

    MatrixXd b_new(n,1);
    ForwardSubSolve forward;
    b_new = forward.solve(D+L_A, b);

    //counter
    int ic = 0;

    //the two vectors in the iteration
    MatrixXd x_old(MatrixXd::Zero(n,1));
    MatrixXd x_new(x0);

    //start the iteration
    //while( (x_new-x_old).norm()>tol && ic <= INT_MAX-1){    //consecutive approx. stopping
    while( r.norm()>stop_iter_resid && ic <= INT_MAX-1){  //residual-based stopping
        x_old = x_new;
        x_new = - forward.solve(D+L_A, U_A*x_old) + b_new;
        //x = - forward_sub(D+L_A, U_A*x) + b_new;
        r = b - A*x_new;
        ic ++;

        //print out the first n iterations
        //if(ic <=3){ cout << "gs " << ic << " iteration: " << endl << x_new << endl << endl; }
    }

    if(ic == INT_MAX-1){
        cout << "gs does not converge" << endl;
    }
    //cout << "gs op count:" << ic << endl;
    return x_new;

}

LinearSolver* GaussSeidel::clone() const{
    return new GaussSeidel(*this);
}

/** SOR **/
Sor::Sor(double w_, double tol_):w(w_), tol(tol_){}

MatrixXd Sor::solve(MatrixXd A, MatrixXd b){
    /*
     * function to solve x for linear equation Ax = b using SOR iteration method
     * converge when A is spd or strictly diagonally dominant
     *
     * input: A: a n*n matrix
     *        b: a n*1 vector
     *        tol: a tolerance factor to stop the iteration
     *        w: choice of a number between 0 and 2 for optimizing the iteration speed
     *
     * output: a n*1 vector as the solution to Ax=b
     *
     */

    //construct the matrix D, L_A, U_A, where L_A + D + U_A = A
    int n = A.rows();

    MatrixXd D(MatrixXd::Zero(n,n));
    MatrixXd L_A(MatrixXd::Zero(n,n));
    MatrixXd U_A(MatrixXd::Zero(n,n));

    for(int i=0; i<n; i++){
        D(i,i) = A(i,i);

        for(int j=0; j<i; j++){
            L_A(i,j) = A(i,j);
        }

        for(int j=i+1; j<n; j++){
            U_A(i,j) = A(i,j);
        }
    }

    //initial guess
    //MatrixXd x0(MatrixXd::Zero(n,1));
    MatrixXd x0(MatrixXd::Constant(n,1,1));

    //stop criteria
    MatrixXd r0(b - A*x0);
    double stop_iter_resid = tol * r0.norm();


    MatrixXd x(n,1);
    x = x0;

    MatrixXd r(n,1);
    r = r0;

    MatrixXd b_new(n,1);
    ForwardSubSolve forward;
    b_new = w * forward.solve(D + w*L_A, b);

    //counter
    int ic = 0;

    //the two vectors in the iteration
    MatrixXd x_old(MatrixXd::Constant(n, 1, 7.77777));
    MatrixXd x_new(x0);

    //start the iteration
    while( (x_new-x_old).norm()>tol && ic <= INT_MAX-1){    //consecutive approx. stopping
    //while( r.norm()>stop_iter_resid && ic <= INT_MAX-1){  //residual-based stopping
        x_old = x_new;
        x_new = forward.solve(D + w*L_A, (1-w)*D*x_old - w*U_A*x_old) + b_new;
        //x = forward_sub(D + w*L_A, (1-w)*D*x - w*U_A*x) + b_new;
        r = b - A*x_new;
        ic ++;

        //print out the first n iterations
        //if(ic <=3){ cout << "sor " << ic << " iteration: " << endl << x_new << endl << endl; }
    }

    if(ic == INT_MAX-1){
        cout << "Jacobi does not converge" << endl;
    }
    //cout << "sor op count:" << ic << endl;
    return x_new;
}

LinearSolver* Sor::clone() const{
    return new Sor(*this);
}

/* Jacobi with Bandwidth m */
JacobiBanded::JacobiBanded(int m_, double tol_):m(m_), tol(tol_){}

MatrixXd JacobiBanded::solve(MatrixXd A, MatrixXd b){
    //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
    //no need for the M, N matrix in this case

    //get the dimension of A
    int n = A.rows();

    //initial guess
    //MatrixXd x0(MatrixXd::Constant(n,1,1));
    MatrixXd x0(MatrixXd::Constant(n,1,0));

    //stop criteria
    MatrixXd r0(b - A*x0);
    double stop_iter_resid = tol * r0.norm();

    MatrixXd r(r0);

    //counter
    int ic = 0;

    //the two vectors in the iteration
    MatrixXd x_old(MatrixXd::Zero(n,1));
    MatrixXd x_new(x0);

    //start the iteration
    //while( (x_new-x_old).norm()>tol && ic <= INT_MAX-1){    //consecutive approx. stopping
    while( r.norm()>stop_iter_resid && ic <= INT_MAX-1){  //residual-based stopping
        x_old = x_new;

        //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
        for(int j=0; j<n; j++){
            //update he entries of x_new one by one
            double x_new_j = 0;
            for(int k=max(j-m+1,0); k<j; k++){ x_new_j = x_new_j + A(j,k)*x_old(k,0); }
            for(int k=j+1; k<min(j+m,n); k++){ x_new_j = x_new_j + A(j,k)*x_old(k,0); }
            x_new(j,0) = (-x_new_j + b(j,0))/A(j,j);
        }

        r = b - A*x_new;
        ic ++;
    }

    if(ic == INT_MAX-1){
        cout << "Jacobi does not converge" << endl;
    }
    cout << endl << "Jacobi op count " << ic << endl;
    return x_new;
}

LinearSolver* JacobiBanded::clone() const{
    return new JacobiBanded(*this);
}


/* Gauss-Seidel with Bandwidth m */
GaussSeidelBanded::GaussSeidelBanded(int m_, double tol_):m(m_), tol(tol_){}

MatrixXd GaussSeidelBanded::solve(MatrixXd A, MatrixXd b){
    //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
    //no need for the M, N matrix in this case

    //get the dimension of A
    int n = A.rows();

    //initial guess
    //MatrixXd x0(MatrixXd::Constant(n,1,1));
    MatrixXd x0(MatrixXd::Constant(n,1,0));

    //stop criteria
    MatrixXd r0(b - A*x0);
    double stop_iter_resid = tol * r0.norm();

    MatrixXd r(r0);

    //counter
    int ic = 0;

    //the two vectors in the iteration
    MatrixXd x_old(MatrixXd::Zero(n,1));
    MatrixXd x_new(x0);

    //start the iteration
    //while( (x_new-x_old).norm()>tol && ic <= INT_MAX-1){    //consecutive approx. stopping
    while( r.norm()>stop_iter_resid && ic <= INT_MAX-1){  //residual-based stopping
        x_old = x_new;

        //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
        for(int j=0; j<n; j++){
            //update the entries of x_new one by one
            double x_new_j = 0;
            for(int k=max(j-m+1,0); k<j; k++){ x_new_j =  x_new_j + A(j,k)*x_new(k,0); }
            for(int k=j+1; k<min(j+m,n); k++){ x_new_j =  x_new_j + A(j,k)*x_old(k,0); }
            x_new(j,0) = (-x_new_j + b(j,0))/A(j,j);
        }

        r = b - A*x_new;
        ic ++;
    }

    if(ic == INT_MAX-1){
        cout << "Gauss-Siedel does not converge" << endl;
    }
    cout << endl << "gs op count " << ic << endl;
    return x_new;

}

LinearSolver* GaussSeidelBanded::clone() const{
    return new GaussSeidelBanded(*this);
}


/* SOR with bandwidth m */
SorBanded::SorBanded(int m_, double w_, double tol_):m(m_), w(w_), tol(tol_){}

MatrixXd SorBanded::solve(MatrixXd A, MatrixXd b){
    //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
    //no need for the M, N matrix in this case

    //get the dimension of A
    int n = A.rows();

    //initial guess
    //MatrixXd x0(MatrixXd::Constant(n,1,1));
    MatrixXd x0(MatrixXd::Constant(n,1,0));

    //stop criteria
    MatrixXd r0(b - A*x0);
    double stop_iter_resid = tol * r0.norm();

    MatrixXd r(r0);

    //counter
    int ic = 0;

    //the two vectors in the iteration
    MatrixXd x_old(MatrixXd::Zero(n,1));
    MatrixXd x_new(x0);

    //start the iteration
    //while( (x_new-x_old).norm()>tol && ic <= INT_MAX-1){    //consecutive approx. stopping
    while( r.norm()>stop_iter_resid && ic <= INT_MAX-1){  //residual-based stopping
        x_old = x_new;

        //since the A is banded, the entry-by-entry version of iteration is more efficient especially when n>>m
        for(int j=0; j<n; j++){
            //update he entries of x_new one by one
            double x_new_j = 0;
            for(int k=max(j-m+1,0); k<j; k++){ x_new_j =  x_new_j + A(j,k)*x_new(k,0); }
            for(int k=j+1; k<min(j+m,n); k++){ x_new_j =  x_new_j + A(j,k)*x_old(k,0); }
            x_new(j,0) = (1-w)*x_old(j,0) + w*( -x_new_j + b(j,0))/A(j,j);
        }

        r = b - A*x_new;
        ic ++;
    }

    if(ic == INT_MAX-1){
        cout << "SOR does not converge" << endl;
    }
    //cout << endl << "sor op count" << ic << endl;
    return x_new;

}

LinearSolver* SorBanded::clone() const{
    return new SorBanded(*this);
}
