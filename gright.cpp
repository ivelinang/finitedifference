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

/* Black Scholes Put */
BlackScholesPutRight::BlackScholesPutRight(double s_, double k_, double vol_, double t_, double r_, double q_):
                                             s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}


double BlackScholesPutRight::get_x_right(){
    double x_right = log(s/k) + (r-q-vol*vol/2.0) * t + 3*vol*sqrt(t);
    return x_right;
}

double BlackScholesPutRight::operator()(double tau){
    return 0;
}

Gright* BlackScholesPutRight::clone() const{
    return new BlackScholesPutRight(*this);
}

/* Black Scholes Call */
BlackScholesCallRight::BlackScholesCallRight(double s_, double k_, double vol_, double t_, double r_, double q_):
                                               s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}


double BlackScholesCallRight::get_x_right(){
    double x_right = log(s/k) + (r-q-vol*vol/2.0) * t + 3*vol*sqrt(t);
    return x_right;
}

double BlackScholesCallRight::operator()(double tau){
    double x_right = get_x_right();

    double a = (r-q) / (vol*vol) - 0.5;
    double b = (r-q) / (vol*vol) + 0.5;
    b = pow(b, 2) + 2*q / (vol*vol);

    double value = k * exp(a*x_right + b*tau) * ( exp(x_right - 2*q*tau / (vol*vol)) - exp(- 2.0*r*tau / (vol*vol)) );
    return value;
}

Gright* BlackScholesCallRight::clone() const{
    return new BlackScholesCallRight(*this);
}
