/*
 * Declaration of the Heat PDE Solver that solves for u(x, tau) in the partial differential equation
 * d(u)/d(tau) = d^2(u)/d(x^2). It will use three methods to solve the PDE:
 *
 *      Forward Euler
 *      Backward Euler
 *      Crank-Nicolson
 *
 * Yongyi Ye
 */

#ifndef HEATPDESOLVER_HPP
#define HEATPDESOLVER_HPP

#include"linearsolver.hpp"
#include"wrapper.hpp"
#include"gleft.hpp"
#include"gright.hpp"
#include"ftau.hpp"

#include<eigen3/Eigen/Dense>
using namespace Eigen;

class HeatPdeSolver{
    /* defining an interface */

    protected:
        double xleft;
        double xright;
        double taufinal;
        Wrapper<Gleft> gleft;
        Wrapper<Gright> gright;
        Wrapper<Ftau> f;

    public:
        /* constructor and destructor */
        HeatPdeSolver(double xleft_, double xright_, double taufinal_,
                        const Gleft &gleft_, const Gright &gright_, const Ftau &f_);
        HeatPdeSolver(const HeatPdeSolver &input);
        virtual ~HeatPdeSolver();
        virtual HeatPdeSolver& operator= (const HeatPdeSolver &input);

        /* the function to solve the pde given boundary conditions by building a (m+1)*(n+1) mesh such that
         * delta-x = (x-right - x-left)/n
         * delta-tau = tau-final/m
         */
        virtual MatrixXd solve_pde(int n, int m) = 0;

        /* virtual copy */
        virtual HeatPdeSolver* clone() const = 0;
};

class ForwardEuler: public HeatPdeSolver{

    public:
        ForwardEuler(double xleft_, double xright_, double taufinal_,
                        const Gleft &gleft_, const Gright &gright_, const Ftau &f_);
        ForwardEuler(const ForwardEuler &input);
        virtual ~ForwardEuler();
        virtual ForwardEuler& operator= (const ForwardEuler &input);

        virtual MatrixXd solve_pde(int n, int m);
        virtual HeatPdeSolver* clone() const;
};

class BackwardEuler: public HeatPdeSolver{

    private:
        Wrapper<LinearSolver> solver;

    public:
        BackwardEuler(double xleft_, double xright_, double taufinal_,
                        const Gleft &gleft_, const Gright &gright_, const Ftau &f_,
                        const LinearSolver &solver_);
        BackwardEuler(const BackwardEuler &input);
        virtual ~BackwardEuler();
        virtual BackwardEuler& operator= (const BackwardEuler &input);

        virtual MatrixXd solve_pde(int n, int m);
        virtual HeatPdeSolver* clone() const;
};

class CrankNicolson: public HeatPdeSolver{

    private:
        Wrapper<LinearSolver> solver;

    public:
        CrankNicolson(double xleft_, double xright_, double taufinal_,
                        const Gleft &gleft_, const Gright &gright_, const Ftau &f_,
                        const LinearSolver &solver_);
        CrankNicolson(const CrankNicolson &input);
        virtual ~CrankNicolson();
        virtual CrankNicolson& operator= (const CrankNicolson &input);

        virtual MatrixXd solve_pde(int n, int m);
        virtual HeatPdeSolver* clone() const;
};




#endif
