
#ifndef MATH9821_HPP
#define MATH9821_HPP

#include<math.h>
using namespace std;

#include<eigen3/Eigen/Dense>

using namespace Eigen;


MatrixXd mult_vector_vector(MatrixXd col_vec, MatrixXd row_vec);
    /*
     * function to compute the column vector - row vector multiplication
     *
     * input: col_vec - column vector 1*m
     *        row_vec - row vector n*1
     *
     * output: the result matrix
     */

MatrixXd mult_matrix_vector(MatrixXd mat, MatrixXd vec);
    /*
     * function to compute the matrix - vector multiplication
     *
     * input: mat - matrix m*n
     *        vec - vector n*1
     *
     * output: the result matrix
     */

void lu_no_pivoting(MatrixXd A, MatrixXd &l, MatrixXd &u);
    /*
     * function to compute the Lu decomposition with NO pivoting of a general n*n matrix such that A = lu
     *
     * input: A - a n*n matrix with all the principal minors det(A(1:i, 1:i)) of A being nonzero
     *        l - reference to lower triangular matrix n*n (will be changed to the l of A = lu)
     *        u - reference to upper triangular matrix n*n (will be changed to the u of A = lu)
     *
     * output: the function returns nothing, but the inputted reference to l, u are changed accordingly
     */


void lu_row_pivoting(MatrixXd &p, MatrixXd A, MatrixXd &l, MatrixXd &u);
    /*
     * function to compute the LU decomposition with row pivoting of a general n*n matrix such that pA = lu
     *
     *input: A - a n*n nonsingular matrix
     *       p - row-permutation matrix to permute the rows of A (a reference to the identity matrix)
     *       l - lower triangular matrix n*n (a reference to n*n zero matrix)
     *       u - upper triangular matrix n*n (a reference to n*n zero matrix)
     *
     *output: the function returns nothing, but the inputted reference to p, l, u are changed accordingly
     *        such that pA = lu
     */



MatrixXd cholesky_general(MatrixXd A);
    /*
     * function to compute the cholesky factor of a (symmetric positive definite) matrix
     *
     * input: A: a n*n symmetric matrix
     *
     * output: a n*n upper triangular matrix as the Cholesky factor of the matrix
     */


MatrixXd cholesky_banded(MatrixXd A, int m);
    /*
     * function to compute the cholesky factor of a (symmetric positive definite) BANDED matrix with bandwidth m
     *
     * input: A: a n*n banded symmetric matrix
     *        m: the bandwidth of the banded matrix A
     *
     * output: a n*n upper triangular matrix as the Cholesky factor of the matrix
     */

#endif

