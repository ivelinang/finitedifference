
/*
 * Implementation for the Boundary Conditions for tau = 0 for the
 * underlying heat pde for solving the Black Sholes pde
 *
 * Yongyi Ye
 */

#include<math.h>
#include<algorithm>
#include"ftau.hpp"
using namespace std;

Ftau::~Ftau(){}

/* Black Scholes Put - Same for European and American */
BsPutTau::BsPutTau(double s_, double k_, double vol_, double t_, double r_, double q_):
                                            s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}

/* get the tau_final for the heat pde such that tau < tau_final */
double BsPutTau::get_tau_final(){
    double tau_final = t * vol*vol / 2.0;
    return tau_final;
}

/* the boundary condition equation */
double BsPutTau::operator()(double x){
    double a = (r-q) / (vol*vol) - 0.5;
    double value = k * exp(a*x) * max<double>(1-exp(x), 0);
    return value;
}

/* virtual copy constructor */
Ftau* BsPutTau::clone() const{
    return new BsPutTau(*this);
}


/* Black Scholes Call - Same for European and American */
BsCallTau::BsCallTau(double s_, double k_, double vol_, double t_, double r_, double q_):
                                            s(s_), k(k_), vol(vol_), t(t_), r(r_), q(q_){}

/* get the tau_final for the heat pde such that tau < tau_final */
double BsCallTau::get_tau_final(){
    double tau_final = t * vol*vol / 2.0;
    return tau_final;
}

/* the boundary condition equation */
double BsCallTau::operator()(double x){
    double a = (r-q) / (vol*vol) - 0.5;
    double value = k * exp(a*x) * max<double>(exp(x)-1, 0);
    return value;
}

/* virtual copy constructor */
Ftau* BsCallTau::clone() const{
    return new BsCallTau(*this);
}
