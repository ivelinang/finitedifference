//implementation of functions for math9821 code library
//
//math 9821 - fall 2012
//yongyi ye
//
//

#include<eigen3/Eigen/Dense>
#include<math.h>
#include<iostream>

#include"decomposer.hpp"

using namespace std;
using namespace Eigen;

//private helper functions

/************************************************************************************************/

int find_i_max(MatrixXd A, int i);
    /*
     * function to find the location of the largest entry (in absolute value) of the vector A(i:n, i)
     * A(i:n, i) is the ith entry of the ith to last row of matrix A
     *
     * input: A: the matrix A
     *        i: starting position A(i,i) where the location of the largest entry is to be found
     *
     * output: the location/index of the largest entry
     */

void switch_row(MatrixXd &mat, int i, int i_max);
    /*
     * function to switch the rows i and i_max of matrix mat
     *
     * input: mat: the reference to n*n matrix whose rows are to be changed
     *        i, i_max: the row numbers of A which are to be interchanged
     *
     * output: returns void, but the reference of mat is changed accordingly
     *
     */
void switch_row_L(MatrixXd &l, int i, int i_max);
    /*
     * function to switch the 1:(i-1) part of the rows i and i_max of l
     *
     * input: l: the reference to n*n lower triangular matrix L in the process of LU decomposition
     *           with row pivoting
     *
     *        i, i_max: the row numbers of l of which parts are to be interchanged
     *
     * output: returns void, but the reference of l is changed accordingly
     *
     */

/*****************************************************************************************************/


MatrixXd mult_vector_vector(MatrixXd col_vec, MatrixXd row_vec){
    //number of rows and columns of the result matrix
    int row = col_vec.rows();
    int col = row_vec.cols();

    //declare the result matrix
    MatrixXd result(row, col);

    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            result(i, j) = col_vec(i,0) * row_vec(0, j);
        }
    }

    return result;
}

MatrixXd mult_matrix_vector(MatrixXd mat, MatrixXd vec){
    //matrix-vector multiplication

    //declare the vector m*1
    MatrixXd result(mat.rows(),1);

    for(int i=0; i<mat.rows(); i++){

        //initialize the entry of the result vector
        result(i,0) = 0;

        for(int j=0; j<mat.cols(); j++){
            result(i,0) = result(i,0) + mat(i,j) * vec(j,0);
        }
    }

    return result;
}

void lu_no_pivoting(MatrixXd A, MatrixXd &l, MatrixXd &u){
    //lu decomposition withOUT row pivoting of a matrix

    //get the number of rows/columns of matrix
    int n = A.rows();

    //initialize l as identity matrices, u as zero matrix
    l = MatrixXd::Identity(n, n);
    u = MatrixXd::Zero(n, n);

    for(int i=0; i<n-1; i++){

        //get the row of u and column of l one by one
        for(int j=i; j<n; j++){
            l(j,i) = A(j,i)/A(i,i);
            u(i,j) = A(i,j);
         }

        //update the matrix A to accommodate recursive calculation for sub_l and sub_u
        for(int j=i+1; j<n; j++){
        for(int k=i+1; k<n; k++){
            A(j,k) = A(j,k) - l(j,i)*u(i,k);
        }
        }

    }

    l(n-1, n-1) = 1; u(n-1, n-1) = A(n-1, n-1);

}


void lu_row_pivoting(MatrixXd &p, MatrixXd A, MatrixXd &l, MatrixXd &u){
    //lu decomposition with row pivoting of a general matrix

    //get the number of rows/columns of matrix
    int n = A.rows();

    //initialize the p and l as identity matrices, u as zero matrix
    p = MatrixXd::Identity(n, n);
    l = MatrixXd::Identity(n, n);
    u = MatrixXd::Zero(n, n);

    for(int i=0; i<n-1; i++){
        //find i_max, the location of the largest entry (in absolute value) of the vector A(i:n, i)
        int i_max = find_i_max(A, i);

        //switch rows i and i_max of A
        switch_row(A, i, i_max);

        //switch the first 1:(i-1) part of rows i and i_max of L
        if(i>=1) { switch_row_L(l, i, i_max); }

        //swith rows i and i_max of p
        switch_row(p, i, i_max);

        for(int j=i; j<n; j++){
            l(j,i) = A(j,i)/A(i,i);
            u(i,j) = A(i,j);
         }

        for(int j=i+1; j<n; j++){
        for(int k=i+1; k<n; k++){
            A(j,k) = A(j,k) - l(j,i)*u(i,k);
        }
        }

    }

    l(n-1, n-1) = 1; u(n-1, n-1) = A(n-1, n-1);

}

/************************* helper functions *************************/

int find_i_max(MatrixXd A, int i){
    /*
     * function to find the location of the largest entry (in absolute value) of the vector A(i:n, i)
     * A(i:n, i) is the ith entry of the ith to last row of matrix A
     *
     * input: A: the matrix A
     *        i: starting position A(i,i) where the location of the largest entry is to be found
     *
     * output: the location/index of the largest entry
     */

    //the number of rows of A
    int row = A.rows();

    //the location of the largest entry
    int i_max = i;

    for(int j=i; j<row; j++){
        if( fabs( A(j,i) )>fabs( A(i_max,i) ) ){
            i_max = j;
        }
    }

    return i_max;
}

void switch_row(MatrixXd &mat, int i, int i_max){
    /*
     * function to switch the rows i and i_max of matrix mat
     *
     * input: mat: the reference to n*n matrix whose rows are to be changed
     *        i, i_max: the row numbers of A which are to be interchanged
     *
     * output: returns void, but the reference of mat is changed accordingly
     *
     */

    int col = mat.cols();

    //temporary row vector 1*col to store the row of mat during the change
    MatrixXd v(1, col);

    //get the contents of row i
    for(int j=0; j<col; j++){ v(0,j) = mat(i,j); }

    //interchange the contents of row i and row i_max
    for(int j=0; j<col; j++){ mat(i,j) = mat(i_max,j); }
    for(int j=0; j<col; j++){ mat(i_max,j) = v(0,j); }

    return;
}

void switch_row_L(MatrixXd &l, int i, int i_max){
    /*
     * function to switch the 1:(i-1) part of the rows i and i_max of l
     *
     * input: l: the reference to n*n lower triangular matrix L in the process of LU decomposition
     *           with row pivoting
     *
     *        i, i_max: the row numbers of l of which parts are to be interchanged
     *
     * output: returns void, but the reference of l is changed accordingly
     *
     */

    //create a temporary vector 1*(i-1) to store part of the matrix l
    MatrixXd ww(1, i);
    for(int j=0; j<=i-1; j++){ ww(0,j) = l(i,j); }

    //assign the content of l(i, 0:(i-1)) to l(i_max, 0:(i-1))
    for(int j=0; j<=i-1; j++){ l(i,j) = l(i_max,j); }

    //assign the content of ww to l(i_max, 0:(i-1))
    for(int j=0; j<=i-1; j++){ l(i_max,j) = ww(0,j); }

    return;
}


/*********************************************************************************/


MatrixXd cholesky_general(MatrixXd A){
    //cholesky decomposition for a symmetric (positive definite) matrix A; A is n*n
    //returns the upper triangular matrix U cholesky factor of A such that U(t)U = A

    //get the dimension of the matrix
    int n = A.rows();

    //the cholesky factor of A
    MatrixXd u(MatrixXd::Zero(n,n));

    for(int i=0; i<n-1; i++){
        u(i,i) = sqrt(A(i,i));

        for(int j=i+1; j<n; j++){
            u(i,j) = A(i,j)/u(i,i);
        }

        for(int j=i+1; j<n; j++){
        for(int k=j; k<n; k++){
            A(j,k) = A(j,k) - u(i,j)*u(i,k);
        }
        }
    }

    u(n-1,n-1) = sqrt( A(n-1,n-1) );

    return u;
}

MatrixXd cholesky_banded(MatrixXd A, int m){
    //cholesky decompsition for a banded n*n symmetric positive definite matrix A with bandwidth m
    //returns the upper triangular matrix U cholesky factor of A such that U(t)U = A

    //get the dimension of the matrix
    int n = A.rows();

    //the cholesky factor of A
    MatrixXd u( MatrixXd::Zero(n,n) );

    for(int i=0; i<n-1; i++){
        u(i,i) = sqrt(A(i,i));

        for(int j=i+1; j<min(i+m,n) ; j++){
            u(i,j) = A(i,j)/u(i,i);
        }

        for(int j=i+1; j<min(i+m,n); j++){
        for(int k=j; k<min(i+m,n); k++){
            A(j,k) = A(j,k) - u(i,j)*u(i,k);
        }
        }
    }

    u(n-1,n-1) = sqrt( A(n-1,n-1) );

    return u;

}

