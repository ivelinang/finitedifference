/*
 * Declaration of the linear solver classes, including:
 *
 *      Forward Substitution
 *      Backward Substitution
 *      LU without pivoting
 *      LU with row pivoting
 *      Cholesky
 *      Jacobi
 *      Gauss-Seidel
 *      SOR
 *
 * they all inherit from the base class interface linearsolver.
 *
 */

#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include<math.h>
using namespace std;

#include<eigen3/Eigen/Dense>
using namespace Eigen;

class LinearSolver{
    public:
        virtual ~LinearSolver();

        // Function to solve x for Ax = b
        virtual MatrixXd solve(MatrixXd A, MatrixXd b) = 0;
        virtual LinearSolver* clone() const = 0;

};

class ForwardSubSolve: public LinearSolver{
    public:
        /*
        * function to compute Ax = b where A is lower triangular
        *
        * input: A - lower triangular matrix n*n
        *        b - column vector 1*n
        *
        * output: the 1*n vector that makes Lx=b
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class BackwardSubSolve: public LinearSolver{
    public:
        /*
        * function to compute Ax = b where A is upper triangular
        *
        * input: A - upper triangular matrix n*n
        *        b - column vector 1*n
        *
        * output: the 1*n vector such that Ux=b
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class LuNoPivSolve: public LinearSolver{
    public:
       /*
        * function to compute x, the solution to the linear system Ax = b without pivoting
        *
        * input: A - a n*n matrix - with nonzero principal minors
        *        b - a n*1 vector
        *
        * output: the n*1 vector x which is the solution to Ax = b
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class LuRowPivSolve: public LinearSolver{
    public:
        /*
        * function to compute x, the solution to the linear system Ax = b with ROW pivoting
        *
        * input: A - a n*n matrix - nonsingular
        *        b - a n*1 vector
        *
        * output: the n*1 vector x which is the solution to Ax = b
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class CholeskySolve: public LinearSolver{
    public:
       /*
        * function to compute x, the solution to the linear system Ax = b with Cholesky factorization
        *
        * input: A - a n*n matrix - spd
        *        b - a n*1 vector
        *
        * output: the n*1 vector x which is the solution to Ax = b
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class CholeskyBandedSolve: public LinearSolver{
    private:
        int bandwidth;
    public:
        CholeskyBandedSolve(int m);
       /*
        * function to compute x, the solution to the linear system Ax = b with Cholesky factorization
        * A is a (symmetric positive definite) BANDED matrix with bandwidth bandwidth
        * input: A - a n*n matrix - spd
        *        b - a n*1 vector
        *        bandwidth: the bandwidth of the banded matrix A
        *
        * output: the n*1 vector x which is the solution to Ax = b
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class Jacobi: public LinearSolver{
    private:
        double tol;
    public:
        Jacobi(double tol_ = pow(10.0, -6.0));
       /*
        * function to solve x for linear equation Ax = b using Jacobi iteration method
        * converge wheren A is strictly diagonally dominant,  weakly diagonally dominant
        * and irreducible, or consistently ordered
        *
        * input: A: a n*n matrix
        *        b: a n*1 vector
        *        tol: a tolerance factor to stop the iteration, default to be 10^-6
        *
        * output: a n*1 vector as the solution to Ax=b
        *
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class GaussSeidel: public LinearSolver{
    private:
        double tol;
    public:
        GaussSeidel(double tol_ = pow(10.0, -6.0));
        /*
         * function to solve x for linear equation Ax = b using gauss-siedel iteration method
         * converge when A is spd, strictly diagonally dominant, diagonally dominant and irreducible,
         * or consistently ordered
         *
         * input: A: a n*n matrix
         *        b: a n*1 vector
         *        tol: a tolerance factor to stop the iteration, default to be 10^-6
         *
         * output: a n*1 vector as the solution to Ax=b
         *
         */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class Sor: public LinearSolver{
    private:
        double tol;
        double w;
    public:
        Sor(double w_, double tol_ = pow(10.0, -6.0));
       /*
        * function to solve x for linear equation Ax = b using SOR iteration method
        * converge when A is spd or strictly diagonally dominant
        *
        * input: A: a n*n matrix
        *        b: a n*1 vector
        *        tol: a tolerance factor to stop the iteration, default to be 10^-6
        *        w: choice of a number between 0 and 2, for optimizing the iteration speed
        *
        *
        * output: a n*1 vector as the solution to Ax=b
        *
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class JacobiBanded: public LinearSolver{
    private:
        int m;
        double tol;
    public:
        JacobiBanded(int m_, double tol_ = pow(10.0, -6.0));
       /*
        * function to solve x for linear equation Ax = b using Jacobi iteration method for banded matrix A with band m
        * converge when A is strictly diagonally dominant,  weakly diagonally dominant and irreducible, or consistently ordered
        *
        * input: A: a n*n matrix
        *        b: a n*1 vector
        *        m: the bandwidth of A
        *        tol: a tolerance factor to stop the iteration, default to be 10^-6
        *
        * output: a n*1 vector as the solution to Ax=b
        *
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class GaussSeidelBanded: public LinearSolver{
    private:
        int m;
        double tol;
    public:
        GaussSeidelBanded(int m_, double tol_ = pow(10.0, -6.0));
       /*
        * function to solve x for linear equation Ax = b using gauss-siedel iteration method for banded matrix A with band m
        * converge when A is spd, strictly diagonally dominant, diagonally dominant and irreducible, or consistently ordered
        *
        * input: A: a n*n matrix
        *        b: a n*1 vector
        *        m: the bandwidth of A
        *        tol: a tolerance factor to stop the iteration, default to be 10^-6
        *
        * output: a n*1 vector as the solution to Ax=b
        *
        */
        virtual MatrixXd solve(MatrixXd A, MatrixXd b);
        virtual LinearSolver* clone() const;
};

class SorBanded: public LinearSolver{
    private:
        int m;
        double w;
        double tol;
    public:
        SorBanded(int m_, double w_ = 1.5, double tol = pow(10.0, -6.0));
    /*
     * function to solve x for linear equation Ax = b using SOR iteration method for banded matrix A with band m
     * converge when A is spd or strictly diagonally dominant
     *
     * input: A: a n*n matrix
     *        b: a n*1 vector
     *        m: the bandwidth of A
     *        tol: a tolerance factor to stop the iteration, default to be 10^-6
     *        w: choice of a number between 0 and 2, default 1.5, for optimizing the iteration speed
     *
     *
     * output: a n*1 vector as the solution to Ax=b
     *
     */
    virtual MatrixXd solve(MatrixXd A, MatrixXd b);
    virtual LinearSolver* clone() const;
};



#endif
