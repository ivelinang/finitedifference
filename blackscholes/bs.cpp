//implementation of the black-sholes class
//
// Yongyi Ye - Jul.29.2012

#include<iostream>
#include<math.h>
#include"bs.hpp"
using namespace std;

/****helper functions*****/

double bs::compute_d1(){
    double d1;
    d1 = log(S/K) + (r - q + sigma*sigma/2.0)*T;
    d1 = d1/(sigma*sqrt(T));
    return d1;
}

double bs::compute_d2(){
    double d2;
    d2 = log(S/K) + (r - q - sigma*sigma/2.0)*T;
    d2 = d2/(sigma*sqrt(T));
    return d2;
}

double bs::normal(double x){
    double z = fabs(x);
    double y = 1/(1+0.2316419*z);
    double a1 = 0.319381530;
    double a2 = -0.356563782;
    double a3 = 1.781477937;
    double a4 = -1.821255978;
    double a5 = 1.330274429;

    double m = 1-exp(-x*x/2.0)*(a1*y+a2*y*y+a3*pow(y,3.0)+a4*pow(y,4.0)+a5*pow(y,5.0))/sqrt(2.0*4*atan(1.0));

    double nn;
    if(x>0){ nn = m;}
    else{ nn = 1-m;}

    return nn;
}

double bs::compute_call_price(){
    double c;
    c = S*exp(-q*T)*normal(d1)-K*exp(-r*T)*normal(d2);
    return c;
}

double bs::compute_put_price(){
    double p;
    p = K*exp(-r*T)*normal(-d2)-S*exp(-q*T)*normal(-d1);
    return p;
}

double bs::compute_call_delta(){
    double delta;
    delta = exp(-q*T)*normal(d1);
    return delta;
}

double bs::compute_put_delta(){
    double delta;
    delta = -exp(-q*T)*normal(-d1);
    return delta;
}

/*****************CHECK****************/
double bs::compute_gamma(){
    double gam;
    gam = exp(-q*T)/( S*sigma*sqrt(T) ) * exp(-d1*d1/2.0) / sqrt(2*4*atan(1));
    return gam;
}
/****************CHECK***************/

double bs::compute_vega(){
    double veg;
    veg = S*exp(-q*T)*sqrt(T)*exp(-d1*d1/2.0)/sqrt(2*4*atan(1));
    return veg;
}

double bs::compute_call_theta(){
  double the;
  the = - S*sigma*exp(-q*T)/(2.0 * sqrt(2.0*4.0*atan(1)*T))*exp(-d1*d1/2.0)
    +q*S*exp(-q*T)*normal(d1)-r*K*exp(-r*T)*normal(d2);
  return the;
}

double bs::compute_put_theta(){
  double the;
  the = - S*sigma*exp(-q*T)/(2.0 * sqrt(2.0*4.0*atan(1)*T))*exp(-d1*d1/2.0)
    -q*S*exp(-q*T)*normal(-d1)+r*K*exp(-r*T)*normal(-d2);
  return the;
}

/*****end*****/


//default constructor
bs::bs():S(0),K(0),T(0),sigma(0),r(0),q(0),d1(0),d2(0),c_price(0),p_price(0),c_delta(0),p_delta(0),gamma_(0),vega_(0){}

//constructor with parameters
bs::bs(double s, double k, double t, double vol, double rate=0, double div=0):S(s), K(k), T(t), sigma(vol), r(rate), q(div){
    //initialize other parameters
    d1 = compute_d1();
    d2 = compute_d2();

    c_price = compute_call_price();
    p_price = compute_put_price();

    c_delta = compute_call_delta();
    p_delta = compute_put_delta();

    gamma_ = compute_gamma();
    vega_ = compute_vega();

    c_theta = compute_call_theta();
    p_theta = compute_put_theta();
}

//getter functions to get the parameters of the option

//call price
double bs::call_price(){
    return c_price;
}

//put price
double bs::put_price(){
    return p_price;
}

//call delta
double bs::call_delta(){
    return c_delta;
}

//put delta
double bs::put_delta(){
    return p_delta;
}

//gamma
double bs::gamma(){
    return gamma_;
}

//vega
double bs::vega(){
    return vega_;
}

//call theta
double bs::call_theta(){
    return c_theta;
}

//put theta
double bs::put_theta(){
    return p_theta;
}
