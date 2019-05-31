/*
 * Implementation of the make functions.
 *
 * Yongyi Ye
 */


#include"../blackscholes/bspdesolver.hpp"
#include"../blackscholes/gleft.hpp"
#include"../blackscholes/gright.hpp"
#include"../blackscholes/ftau.hpp"
#include"../blackscholes/checkearlyexercise.hpp"

#include"../heatpdesolvers/heatpdesolver.hpp"
#include"../linearalgebra/linearsolver.hpp"
#include"makesolvers.hpp"


/* European Call */
BsPdeSolver make_euro_call_fe(double s, double k, double vol, double t, double r, double q){
    BsEuropeanCallLeft left(s, k, vol, t, r, q);
    BsEuropeanCallRight right(s, k, vol, t, r, q);
    BsCallTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    ForwardEuler fe(xleft, xright, taufinal, left, right, f);
    BsPdeSolver result(s, k, vol, t, r, q, fe);
    return result;
}

BsPdeSolver make_euro_call_be_lu(double s, double k, double vol, double t, double r, double q) {
	BsEuropeanCallLeft left(s, k, vol, t, r, q);
	BsEuropeanCallRight right(s, k, vol, t, r, q);
	BsCallTau f(s, k, vol, t, r, q);

	double xleft = left.get_x_left();
	double xright = right.get_x_right();
	double taufinal = f.get_tau_final();

	LuNoPivSolve LuSolve ;
	BackwardEuler be(xleft, xright, taufinal, left, right, f, LuSolve);
	BsPdeSolver result(s, k, vol, t, r, q, be);
	return result;
}

BsPdeSolver make_euro_call_be_cholesky(double s, double k, double vol, double t, double r, double q){
    BsEuropeanCallLeft left(s, k, vol, t, r, q);
    BsEuropeanCallRight right(s, k, vol, t, r, q);
    BsCallTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    CholeskySolve cholesky;
    BackwardEuler be(xleft, xright, taufinal, left, right, f, cholesky);
    BsPdeSolver result(s, k, vol, t, r, q, be);
    return result;
}

BsPdeSolver make_euro_call_be_sor(double s, double k, double vol, double t, double r, double q,
                                    double w, double tol){

    BsEuropeanCallLeft left(s, k, vol, t, r, q);
    BsEuropeanCallRight right(s, k, vol, t, r, q);
    BsCallTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    Sor sorsolver(w, tol);
    BackwardEuler be(xleft, xright, taufinal, left, right, f, sorsolver);
    BsPdeSolver result(s, k, vol, t, r, q, be);
    return result;
}

BsPdeSolver make_euro_call_cn_cholesky(double s, double k, double vol, double t, double r, double q){
    BsEuropeanCallLeft left(s, k, vol, t, r, q);
    BsEuropeanCallRight right(s, k, vol, t, r, q);
    BsCallTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    CholeskySolve cholesky;
    CrankNicolson cn(xleft, xright, taufinal, left, right, f, cholesky);
    BsPdeSolver result(s, k, vol, t, r, q, cn);
    return result;
}

BsPdeSolver make_euro_call_cn_sor(double s, double k, double vol, double t, double r, double q,
                                    double w, double tol){

    BsEuropeanCallLeft left(s, k, vol, t, r, q);
    BsEuropeanCallRight right(s, k, vol, t, r, q);
    BsCallTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    Sor sorsolver(w, tol);
    CrankNicolson cn(xleft, xright, taufinal, left, right, f, sorsolver);
    BsPdeSolver result(s, k, vol, t, r, q, cn);
    return result;
}


/* European Put */
BsPdeSolver make_euro_put_fe(double s, double k, double vol, double t, double r, double q){
    BsEuropeanPutLeft left(s, k, vol, t, r, q);
    BsEuropeanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    ForwardEuler fe(xleft, xright, taufinal, left, right, f);
    BsPdeSolver result(s, k, vol, t, r, q, fe);
    return result;
}

BsPdeSolver make_euro_put_be_cholesky(double s, double k, double vol, double t, double r, double q){
    BsEuropeanPutLeft left(s, k, vol, t, r, q);
    BsEuropeanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    CholeskySolve cholesky;
    BackwardEuler be(xleft, xright, taufinal, left, right, f, cholesky);
    BsPdeSolver result(s, k, vol, t, r, q, be);
    return result;
}

BsPdeSolver make_euro_put_be_sor(double s, double k, double vol, double t, double r, double q,
                                    double w, double tol){

    BsEuropeanPutLeft left(s, k, vol, t, r, q);
    BsEuropeanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    Sor sorsolver(w, tol);
    BackwardEuler be(xleft, xright, taufinal, left, right, f, sorsolver);
    BsPdeSolver result(s, k, vol, t, r, q, be);
    return result;
}

BsPdeSolver make_euro_put_cn_cholesky(double s, double k, double vol, double t, double r, double q){
    BsEuropeanPutLeft left(s, k, vol, t, r, q);
    BsEuropeanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    CholeskySolve cholesky;
    CrankNicolson cn(xleft, xright, taufinal, left, right, f, cholesky);
    BsPdeSolver result(s, k, vol, t, r, q, cn);
    return result;
}

BsPdeSolver make_euro_put_cn_sor(double s, double k, double vol, double t, double r, double q,
                                    double w, double tol){

    BsEuropeanPutLeft left(s, k, vol, t, r, q);
    BsEuropeanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    Sor sorsolver(w, tol);
    CrankNicolson cn(xleft, xright, taufinal, left, right, f, sorsolver);
    BsPdeSolver result(s, k, vol, t, r, q, cn);
    return result;
}


/* American Call */
BsPdeSolver make_amer_call_fe(double s, double k, double vol, double t, double r, double q){
    BsAmericanCallLeft left(s, k, vol, t, r, q);
    BsAmericanCallRight right(s, k, vol, t, r, q);
    BsCallTau f(s, k, vol, t, r, q);
    CheckAmericanCall checker(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    EarlyExForwardEuler exfe(xleft, xright, taufinal, left, right, f, checker);
    BsPdeSolver result(s, k, vol, t, r, q, exfe);
    return result;
}

BsPdeSolver make_amer_call_cn(double s, double k, double vol, double t, double r, double q,
                                double w, double tol){

    BsAmericanCallLeft left(s, k, vol, t, r, q);
    BsAmericanCallRight right(s, k, vol, t, r, q);
    BsCallTau f(s, k, vol, t, r, q);
    CheckAmericanCall checker(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    EarlyExCrankNicolson excn(xleft, xright, taufinal, left, right, f, checker, w, tol);
    BsPdeSolver result(s, k, vol, t, r, q, excn);
    return result;
}


/* American Put */
BsPdeSolver make_amer_put_fe(double s, double k, double vol, double t, double r, double q){
    BsAmericanPutLeft left(s, k, vol, t, r, q);
    BsAmericanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);
    CheckAmericanPut checker(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    EarlyExForwardEuler exfe(xleft, xright, taufinal, left, right, f, checker);
    BsPdeSolver result(s, k, vol, t, r, q, exfe);
    return result;
}

BsPdeSolver make_amer_put_cn(double s, double k, double vol, double t, double r, double q,
                                double w, double tol){

    BsAmericanPutLeft left(s, k, vol, t, r, q);
    BsAmericanPutRight right(s, k, vol, t, r, q);
    BsPutTau f(s, k, vol, t, r, q);
    CheckAmericanPut checker(s, k, vol, t, r, q);

    double xleft = left.get_x_left();
    double xright = right.get_x_right();
    double taufinal = f.get_tau_final();

    EarlyExCrankNicolson excn(xleft, xright, taufinal, left, right, f, checker, w, tol);
    BsPdeSolver result(s, k, vol, t, r, q, excn);
    return result;
}
