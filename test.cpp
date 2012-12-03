/*
 * Test and homeword file for the finite difference partial differential equation
 * solver classes.
 *
 */


#include<iostream>
#include<eigen3/Eigen/Dense>

#include"linearsolver.hpp"
#include"heatpdesolver.hpp"
#include"gleft.hpp"
#include"gright.hpp"
#include"ftau.hpp"

using namespace std;
using namespace Eigen;

int main(){
    cout.precision(9);

    /* the boundary conditions */
    Hw8left left;
    Hw8right right;
    Hw8tau f;

    /* solvers to be used */
    CholeskySolve cholesky;
    Sor sor12(1.2);

    /* the Heat PDE solvers */
    ForwardEuler frd(-2.0, 2.0, 1.0, left, right, f);

    BackwardEuler bdcho(-2.0, 2.0, 1.0, left, right, f, cholesky);
    BackwardEuler bdsor(-2.0, 2.0, 1.0, left, right, f, sor12);

    CrankNicolson cncho(-2.0, 2.0, 1.0, left, right, f, cholesky);
    CrankNicolson cnsor(-2.0, 2.0, 1.0, left, right, f, sor12);

    /* solve */

    // Forward Euler
    cout << "forward euler alpha 0.125" << endl;
    cout << frd.solve_pde(4, 8) << endl << endl;

    cout << "forward euler alpha 0.5" << endl;
    cout << frd.solve_pde(8, 8) << endl << endl;

    cout << "forward euler alpha 2" << endl;
    cout << frd.solve_pde(16, 8) << endl << endl;


    // Backward Euler + Cholesky
    cout << "backward euler alpha 0.125 + cholesky" << endl;
    cout << bdcho.solve_pde(4, 8) << endl << endl;

    cout << "backward euler alpha 0.5 + cholesky" << endl;
    cout << bdcho.solve_pde(8, 8) << endl << endl;

    cout << "backward euler alpha 2 + cholesky" << endl;
    cout << bdcho.solve_pde(16, 8) << endl << endl;


    // Backward Euler + SOR (w = 1.2)
    cout << "backward euler alpha 0.125 + sor" << endl;
    cout << bdsor.solve_pde(4, 8) << endl << endl;

    cout << "backward euler alpha 0.5 + sor" << endl;
    cout << bdsor.solve_pde(8, 8) << endl << endl;

    cout << "backward euler alpha 2 + sor" << endl;
    cout << bdsor.solve_pde(16, 8) << endl << endl;

    //Crank-Nicolson + Cholesky
    cout << "crank-nicolson alpha 0.125 + cholesky" << endl;
    cout << cncho.solve_pde(4, 8) << endl << endl;

    cout << "crank-nicolson alpha 0.5 + cholesky" << endl;
    cout << cncho.solve_pde(8, 8) << endl << endl;

    cout << "crank-nicolson alpha 2 + cholesky" << endl;
    cout << cncho.solve_pde(16, 8) << endl << endl;


    //Crank-Nicolson + SOR (w = 1.2)
    cout << "crank-nicolson alpha 0.125 + sor" << endl;
    cout << cnsor.solve_pde(4, 8) << endl << endl;

    cout << "crank-nicolson alpha 0.5 + sor" << endl;
    cout << cnsor.solve_pde(8, 8) << endl << endl;

    cout << "crank-nicolson alpha 2 + sor" << endl;
    cout << cnsor.solve_pde(16, 8) << endl << endl;



}
