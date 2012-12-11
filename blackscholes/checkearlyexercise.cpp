/*
 * Implementation of the check-early-exercise classes
 *
 * Yongyi Ye
 */

#include<math.h>
#include<algorithm>
#include"checkearlyexercise.hpp"

using namespace std;

CheckEarlyExercise::CheckEarlyExercise(double s_, double k_, double vol_, double t_, double r_, double q_):
                                        s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){

    // x_domain = x_right - x_left = 6 * vol * sqrt(t)
    x_left = log(s/k) + (r - q - vol*vol/2.0) * t - 3.0 * vol * sqrt(t);
    x_domain = 6 * vol * sqrt(t);
    // tau_domain = tau_final - 0 = tau_final
    tau_domain = t * vol * vol / 2.0;

    // the a and b involved in the transformation from black-scholes pde to heat pde
    a = (r-q) / (vol*vol) - 0.5;
    b = pow((r-q) / (vol*vol) + 0.5, 2) + (2.0 * q) / (vol*vol);
}

CheckEarlyExercise::CheckEarlyExercise(const CheckEarlyExercise &input):s(input.s), k(input.k), vol(input.vol), t(input.t),
                                                                        r(input.r), q(input.q), x_left(input.x_left),
                                                                        x_domain(input.x_domain), tau_domain(input.tau_domain),
                                                                        a(input.a), b(input.b){}

CheckEarlyExercise::~CheckEarlyExercise(){}

CheckEarlyExercise& CheckEarlyExercise::operator= (const CheckEarlyExercise &input){
    s = input.s;
    k = input.k;
    vol = input.vol;
    t = input.t;
    r = input.r;
    q = input.q;
    x_left = input.x_left;
    x_domain = input.x_domain;
    tau_domain = input.tau_domain;
    a = input.a;
    b = input.b;

    return *this;
}


/* check early exercise for American put */
CheckAmericanPut::CheckAmericanPut(double s_, double k_, double vol_, double t_, double r_, double q_):
                                      CheckEarlyExercise(s_, k_, vol_, t_, r_, q_){}

CheckAmericanPut::CheckAmericanPut(const CheckAmericanPut &input): CheckEarlyExercise(input){}

CheckAmericanPut::~CheckAmericanPut(){}

CheckAmericanPut& CheckAmericanPut::operator= (const CheckAmericanPut &input){
    CheckEarlyExercise::operator=(input);
    return *this;
}

double CheckAmericanPut::check_early_exercise(int N, int M, int n, int m){
    // find the value pair at point (m, n)
    double delta_x = x_domain/static_cast<double>(N);
    double x_n = x_left + n * delta_x;

    double delta_tau = tau_domain/static_cast<double>(M);
    double tau_m = m * delta_tau;

    double check = k * exp(a*x_n + b*tau_m) * max<double>(1-exp(x_n), 0);
    return check;
}

CheckEarlyExercise* CheckAmericanPut::clone() const{
    return new CheckAmericanPut(*this);
}

/* check early exercise for American call */
CheckAmericanCall::CheckAmericanCall(double s_, double k_, double vol_, double t_, double r_, double q_):
                                        CheckEarlyExercise(s_, k_, vol_, t_, r_, q_){}

CheckAmericanCall::CheckAmericanCall(const CheckAmericanCall &input):CheckEarlyExercise(input){}

CheckAmericanCall::~CheckAmericanCall(){}

CheckAmericanCall& CheckAmericanCall::operator= (const CheckAmericanCall &input){
    CheckEarlyExercise::operator=(input);
    return *this;
}

double CheckAmericanCall::check_early_exercise(int N, int M, int n, int m){
    // find the value pair at point (m, n)
    double delta_x = x_domain/static_cast<double>(N);
    double x_n = x_left + n * delta_x;

    double delta_tau = tau_domain/static_cast<double>(M);
    double tau_m = m * delta_tau;

    double check = k * exp(a*x_n + b*tau_m) * max<double>(exp(x_n)-1, 0);
    return check;
}

CheckEarlyExercise* CheckAmericanCall::clone() const{
    return new CheckAmericanCall(*this);
}
