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
#include"bspdesolver.hpp"

using namespace std;
using namespace Eigen;

int main(){
    cout.precision(11);

    /* Code Test 4 */
    double s = 40;
    double k = 42;
    double vol = 0.3;
    double t = 0.5;
    double r = 0.03;
    double q = 0.01;

    BsEuropeanPutLeft left(s, k, vol, t, r, q);
    BsEuropeanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k , vol, t, r, q);
    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double tau = f.get_tau_final();
    ForwardEuler frd(xleft, xright, tau, left, right, f);
    BsPdeSolver efe(s, k, vol, t, r, q, frd);

    BsAmericanPutLeft aleft(s, k, vol, t, r, q);
    BsAmericanPutRight aright(s, k, vol, t, r, q);
    CheckAmericanPut checker(s, k, vol, t, r, q);
    double axleft = aleft.get_x_left();
    double axright = aright.get_x_right();
    double atau = f.get_tau_final();
    EarlyExForwardEuler putfe(axleft, axright, atau, aleft, aright, f, checker);
    BsPdeSolver afe(s, k, vol, t, r, q, putfe);

    double al = 0.5;
    int m = 4;

    cout << "N " << efe.compute_n(al, m) << endl
         << "xleft " << xleft << " " << axleft << endl
         << "xright" << xright << " " << axright << endl << endl;

    cout << "European:" << endl;
    cout << efe.solve_pde(al, m) << endl << endl;
    cout << "price: " << efe.compute_price(al, m) << endl;
    cout << "delta " << efe.compute_delta(al, m) << endl;
    cout << "gamma " << efe.compute_gamma(al, m) << endl << endl;

    cout << "American:" << endl;
    cout << afe.solve_pde(al, m) << endl << endl;
    cout << "price: " << afe.compute_price(al, m) << endl;
    cout << "delta " << afe.compute_delta(al, m) << endl;
    cout << "gamma " << afe.compute_gamma(al, m) << endl;


    /* Homework 10 */
    /*
    double s = 41;
    double k = 40;
    double vol = 0.35;
    double t = 0.75;
    double r = 0.04;
    double q = 0.02;

    BsEuropeanPutLeft left(s, k, vol, t, r, q);
    BsEuropeanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double tau = f.get_tau_final();

    CholeskySolve cholesky;
    Sor sor12(1.2);

    ForwardEuler frd(xleft, xright, tau, left, right, f);
    BackwardEuler bdcho(xleft, xright, tau, left, right, f, cholesky);
    CrankNicolson cncho(xleft, xright, tau, left, right, f, cholesky);
    BsPdeSolver bfe(s, k, vol, t, r, q, frd);
    BsPdeSolver bbe(s, k, vol, t, r, q, bdcho);
    BsPdeSolver bcn(s, k, vol, t, r, q, cncho);
    cout << "Forward Euler European Put" << endl << bfe.compute_price(0.45, 4) << endl << endl;
    cout << "Forward Euler" << endl << bfe.solve_pde(0.45, 4) << endl << endl;
    cout << "Backward Euler" << endl << bbe.solve_pde(0.45, 4) << endl << endl;
    cout << "Crank Nicolson" << endl << bcn.solve_pde(0.45, 4) << endl << endl;

    BsAmericanPutLeft aleft(s, k, vol, t, r, q);
    BsAmericanPutRight aright(s, k, vol, t, r, q);
    CheckAmericanPut checker(s, k, vol, t, r, q);
    double axleft = aleft.get_x_left();
    double axright = aright.get_x_right();
    double atau = f.get_tau_final();

    EarlyExForwardEuler aputfe(axleft, axright, atau, aleft, aright, f, checker);
    EarlyExCrankNicolson aputcn(axleft, axright, atau, aleft, aright, f, checker);
    BsPdeSolver afe(s, k, vol, t, r, q, aputfe);
    BsPdeSolver acn(s, k, vol, t, r, q, aputcn);

    cout << "American Put FE" << endl << afe.solve_pde(0.45, 4) << endl << "price" << endl << afe.compute_price(0.45, 4) << endl << endl;
    cout << "American Put CN" << endl << acn.solve_pde(0.45, 4) << endl << "price" << endl << acn.compute_price(0.45, 4) << endl << endl;

    */


    /* Homework 8 */
    /*************
    // the boundary conditions
    Hw8left left;
    Hw8right right;
    Hw8tau f;

    // solvers to be used
    CholeskySolve cholesky;
    Sor sor12(1.2);

    // the Heat PDE solvers
    ForwardEuler frd(-2.0, 2.0, 1.0, left, right, f);
    BackwardEuler bdcho(-2.0, 2.0, 1.0, left, right, f, cholesky);
    BackwardEuler bdsor(-2.0, 2.0, 1.0, left, right, f, sor12);
    CrankNicolson cncho(-2.0, 2.0, 1.0, left, right, f, cholesky);
    CrankNicolson cnsor(-2.0, 2.0, 1.0, left, right, f, sor12);

    // solve
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
    *********/


}
