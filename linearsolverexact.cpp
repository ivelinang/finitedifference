/*
 * Implementation of the linear solver classes that uses decomposition to get
 * solutions, including:
 *
 *      Forward Substitution
 *      Backward Substitution
 *      LU without pivoting
 *      LU with row pivoting
 *      Cholesky
 *
 * they all inherit from the base class interface linearsolver.
 *
 */

#include<math.h>
#include<eigen3/Eigen/Dense>

#include"linearsolver.hpp"
#include"decomposer.hpp"

using namespace std;
using namespace Eigen;

/* Linear Solver */
LinearSolver::~LinearSolver(){}


/* Forward Sub Solve */
MatrixXd ForwardSubSolve::solve(MatrixXd A, MatrixXd b){
    //forward substitution
    int n = A.rows();

    //declare the x vector
    MatrixXd x(n,1);
    x(0,0) = b(0,0)/A(0,0);

    for(int i=1; i<n; i++){
        double sum = 0;

        for(int j=0; j <= i-1; j++){
            sum = sum + A(i,j)*x(j,0);
        }
        x(i) = ( b(i,0) - sum)/A(i,i);
    }
    return x;
}

LinearSolver* ForwardSubSolve::clone() const{
    return new ForwardSubSolve(*this);
}


/* Backward Sub Solve */
MatrixXd BackwardSubSolve::solve(MatrixXd A, MatrixXd b){
    //backward substitution
    int n = A.rows();

    //declare the x vector
    MatrixXd x(n,1);
    x(n-1, 0) = b(n-1,0)/A(n-1,n-1);

    for(int i=n-2; i>=0; i--){
        double sum = 0;

        for(int j=i+1; j<n; j++){
            sum = sum + A(i,j)*x(j,0);
        }
        x(i,0) = ( b(i,0) -sum ) /A(i,i);
    }

    return x;
}

LinearSolver* BackwardSubSolve::clone() const{
    return new BackwardSubSolve(*this);
}


/* LU Decomposition with NO pivoting */
MatrixXd LuNoPivSolve::solve(MatrixXd A, MatrixXd b){
    // linear solve function to solve x for the linear system Ax=b
    // by using lu decomposition with no pivoting

    //get the dimension of A
    int n = A.rows();

    //create to zero n*n matrices l and u
    MatrixXd l( MatrixXd::Zero(n,n) );
    MatrixXd u( MatrixXd::Zero(n,n) );

    lu_no_pivoting(A, l, u);

    ForwardSubSolve forward;
    BackwardSubSolve backward;

    MatrixXd y( forward.solve(l,b) );
    MatrixXd x( backward.solve(u,y) );

    return x;
}

LinearSolver* LuNoPivSolve::clone() const{
    return new LuNoPivSolve(*this);
}


/* LU Decomposition with row pivoting */
MatrixXd LuRowPivSolve::solve(MatrixXd A, MatrixXd b){
    //linear solve function to solve x for the linear system Ax=b by using lu decomposition with row pivoting

    //get the dimension of A
    int n = A.rows();

    //create identity/zero n*n matrices p, l and u
    MatrixXd p( MatrixXd::Identity(n,n) );
    MatrixXd l( MatrixXd::Zero(n,n) );
    MatrixXd u( MatrixXd::Zero(n,n) );

    lu_row_pivoting(p, A, l, u);
    ForwardSubSolve forward;
    BackwardSubSolve backward;

    MatrixXd y( forward.solve(l,p*b) );
    MatrixXd x( backward.solve(u,y) );

    return x;
}

LinearSolver* LuRowPivSolve::clone() const{
    return new LuRowPivSolve(*this);
}


/* Linear Solve with Cholesky Decomposition */
MatrixXd CholeskySolve::solve(MatrixXd A, MatrixXd b){
    // solve Ax = b using Cholesky decomposition

    //get the dimension of A
    int n = A.rows();

    //run cholesky to get the cholesky factor
    MatrixXd u(cholesky_general(A));

    ForwardSubSolve forward;
    BackwardSubSolve backward;

    MatrixXd y( forward.solve(u.transpose(),b) );
    MatrixXd x( backward.solve(u,y) );

    return x;
}

LinearSolver* CholeskySolve::clone() const{
    return new CholeskySolve(*this);
}

/* Linear Solve with Cholesky Decomposition of banded Matrix A with bandwith m */
CholeskyBandedSolve::CholeskyBandedSolve(int m): bandwidth(m){}

MatrixXd CholeskyBandedSolve::solve(MatrixXd A, MatrixXd b){
    // solve Ax = b using Cholesky decomposition

    //get the dimension of A
    int n = A.rows();

    //run cholesky to get the cholesky factor
    MatrixXd u(cholesky_banded(A, bandwidth));

    ForwardSubSolve forward;
    BackwardSubSolve backward;

    MatrixXd y( forward.solve(u.transpose(),b) );
    MatrixXd x( backward.solve(u,y) );

    return x;
}

LinearSolver* CholeskyBandedSolve::clone() const{
    return new CholeskyBandedSolve(*this);
}
