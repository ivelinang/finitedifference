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

/* Black Scholes Put */
BlackScholesPutLeft::BlackScholesPutLeft(double s_, double k_, double vol_, double t_, double r_, double q_):
                                            s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}


double BlackScholesPutLeft::get_x_left(){
    double x_left = log(s/k) + (r-q-vol*vol/2.0) * t - 3*vol*sqrt(t);
    return x_left;
}

double BlackScholesPutLeft::operator()(double tau){
    double x_left = get_x_left();

    double a = (r-q) / (vol*vol) - 0.5;
    double b = (r-q) / (vol*vol) + 0.5;
    b = pow(b, 2) + 2*q / (vol*vol);

    double value = k * exp(a*x_left + b*tau) * ( exp(- 2.0*r*tau / (vol*vol)) - exp(x_left - 2*q*tau / (vol*vol)) );
    return value;
}

Gleft* BlackScholesPutLeft::clone() const{
    return new BlackScholesPutLeft(*this);
}

/* Black Scholes Call */
BlackScholesCallLeft::BlackScholesCallLeft(double s_, double k_, double vol_, double t_, double r_, double q_):
                                    s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}


double BlackScholesCallLeft::get_x_left(){
    double x_left = log(s/k) + (r-q-vol*vol/2.0) * t - 3*vol*sqrt(t);
    return x_left;
}

double BlackScholesCallLeft::operator()(double tau){
    return 0;
}

Gleft* BlackScholesCallLeft::clone() const{
    return new BlackScholesCallLeft(*this);
}
