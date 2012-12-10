/*
 * Implementation for the right Boundary Conditions for x_right and g_right for the
 * underlying heat pde for solving the Black Sholes pde
 *
 * Yongyi Ye
 */

#include<math.h>
#include"gright.hpp"
using namespace std;

Gright::~Gright(){}

/* Black Scholes European Put */
BsEuropeanPutRight::BsEuropeanPutRight(double s_, double k_, double vol_, double t_, double r_, double q_):
                                             s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}


double BsEuropeanPutRight::get_x_right(){
    double x_right = log(s/k) + (r-q-vol*vol/2.0) * t + 3*vol*sqrt(t);
    return x_right;
}

double BsEuropeanPutRight::operator()(double tau){
    return 0;
}

Gright* BsEuropeanPutRight::clone() const{
    return new BsEuropeanPutRight(*this);
}


/* Black Scholes European Call */
BsEuropeanCallRight::BsEuropeanCallRight(double s_, double k_, double vol_, double t_, double r_, double q_):
                                               s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){

    a = (r-q) / (vol*vol) - 0.5;
    b = (r-q) / (vol*vol) + 0.5;
    b = pow(b, 2) + 2*q / (vol*vol);
}


double BsEuropeanCallRight::get_x_right(){
    double x_right = log(s/k) + (r-q-vol*vol/2.0) * t + 3*vol*sqrt(t);
    return x_right;
}

double BsEuropeanCallRight::operator()(double tau){
    double x_right = get_x_right();

    double value = k * exp(a*x_right + b*tau) * ( exp(x_right - 2*q*tau / (vol*vol)) - exp(- 2.0*r*tau / (vol*vol)) );
    return value;
}

Gright* BsEuropeanCallRight::clone() const{
    return new BsEuropeanCallRight(*this);
}


/* Black Scholes American Put */
BsAmericanPutRight::BsAmericanPutRight(double s_, double k_, double vol_, double t_, double r_, double q_):
                                             s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}


double BsAmericanPutRight::get_x_right(){
    // same as European put
    double x_right = log(s/k) + (r-q-vol*vol/2.0) * t + 3*vol*sqrt(t);
    return x_right;
}

double BsAmericanPutRight::operator()(double tau){
    // as s-> inf, option has 0 value
    return 0;
}

Gright* BsAmericanPutRight::clone() const{
    return new BsAmericanPutRight(*this);
}

/* Black Scholes American Call */
BsAmericanCallRight::BsAmericanCallRight(double s_, double k_, double vol_, double t_, double r_, double q_):
                                               s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){

    a = (r-q) / (vol*vol) - 0.5;
    b = (r-q) / (vol*vol) + 0.5;
    b = pow(b, 2) + 2*q / (vol*vol);
}


double BsAmericanCallRight::get_x_right(){
    // same as European call
    double x_right = log(s/k) + (r-q-vol*vol/2.0) * t + 3*vol*sqrt(t);
    return x_right;
}

double BsAmericanCallRight::operator()(double tau){
    double x_right = get_x_right();

    // as s->inf, it is optimal to exercise early, therefore boundary condition is s-k, or as follows
    // in terms of x_right and tau
    double value = k * exp(a*x_right + b*tau) * (exp(x_right) - 1);
    return value;
}

Gright* BsAmericanCallRight::clone() const{
    return new BsAmericanCallRight(*this);
}

