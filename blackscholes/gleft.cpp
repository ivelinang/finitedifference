/*
 * Implementation for the left Boundary Conditions for x_left and g_left for the
 * underlying heat pde for solving the Black Sholes pde
 *
 * Yongyi Ye
 */

#include<math.h>
#include"gleft.hpp"
using namespace std;

Gleft::~Gleft(){}

/* Black Scholes European Put */
BsEuropeanPutLeft::BsEuropeanPutLeft(double s_, double k_, double vol_, double t_, double r_, double q_):
                                            s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){
    a = (r-q) / (vol*vol) - 0.5;
    b = (r-q) / (vol*vol) + 0.5;
    b = pow(b, 2) + 2*q / (vol*vol);
}


double BsEuropeanPutLeft::get_x_left(){
    double x_left = log(s/k) + (r-q-vol*vol/2.0) * t - 3*vol*sqrt(t);
    return x_left;
}

double BsEuropeanPutLeft::operator()(double tau){
    double x_left = get_x_left();

    double value = k * exp(a*x_left + b*tau) * ( exp(- 2.0*r*tau / (vol*vol)) - exp(x_left - 2*q*tau / (vol*vol)) );
    return value;
}

Gleft* BsEuropeanPutLeft::clone() const{
    return new BsEuropeanPutLeft(*this);
}


/* Black Scholes European Call */
BsEuropeanCallLeft::BsEuropeanCallLeft(double s_, double k_, double vol_, double t_, double r_, double q_):
                                    s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}



double BsEuropeanCallLeft::get_x_left(){
    double x_left = log(s/k) + (r-q-vol*vol/2.0) * t - 3*vol*sqrt(t);
    return x_left;
}

double BsEuropeanCallLeft::operator()(double tau){
    return 0;
}

Gleft* BsEuropeanCallLeft::clone() const{
    return new BsEuropeanCallLeft(*this);
}


/* Black Scholes American Put */
BsAmericanPutLeft::BsAmericanPutLeft(double s_, double k_, double vol_, double t_, double r_, double q_):
                                            s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){
    a = (r-q) / (vol*vol) - 0.5;
    b = (r-q) / (vol*vol) + 0.5;
    b = pow(b, 2) + 2*q / (vol*vol);
}


double BsAmericanPutLeft::get_x_left(){
    // same as European put
    double x_left = log(s/k) + (r-q-vol*vol/2.0) * t - 3*vol*sqrt(t);
    return x_left;
}

double BsAmericanPutLeft::operator()(double tau){
    double x_left = get_x_left();

    // as s->0, it is optimal to exercise early, therefore boundary condition is k-s, or as follows
    // in terms of x_left and tau
    double value = k * exp(a*x_left + b*tau) * (1 - exp(x_left));
    return value;
}

Gleft* BsAmericanPutLeft::clone() const{
    return new BsAmericanPutLeft(*this);
}

/* Black Scholes American Call */
BsAmericanCallLeft::BsAmericanCallLeft(double s_, double k_, double vol_, double t_, double r_, double q_):
                                    s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}



double BsAmericanCallLeft::get_x_left(){
    // same as European call
    double x_left = log(s/k) + (r-q-vol*vol/2.0) * t - 3*vol*sqrt(t);
    return x_left;
}

double BsAmericanCallLeft::operator()(double tau){
    // as s->0, option has 0 value
    return 0;
}

Gleft* BsAmericanCallLeft::clone() const{
    return new BsAmericanCallLeft(*this);
}



