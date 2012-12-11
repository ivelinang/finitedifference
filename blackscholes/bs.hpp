//declaration of the Black-Scholes equation class
//input: the spot price, strike, maturity (in years), volatility, risk-free
//	interest rate, and divident yield
//it has all the getter functions to get the values of the price and greeks for
//	the option (either call or put)
//
// Yongyi Ye - Jul.29.2012

#ifndef BS_HPP
#define BS_HPP

class bs{

private:
    double S;   //spot price
    double K;   //strike
    double T;   //maturity in years
    double sigma;   //volatility
    double r;   //riskfree interest rate
    double q;   //continuous divident rate

    double d1;  //the B-S d1 for computation in standard normal
    double d2;  //the B-S d2 for computation in standard normal

    double c_price; //call price
    double p_price; //put price

    double c_delta; //call delta
    double p_delta; //put delta

    double gamma_;   //gamma for both call and put
    double vega_;    //vega for both call and put

    double c_theta; //call theta
    double p_theta; //put theta

    /***helper functions***/

    double compute_d1();//compute the value for d1
    double compute_d2();//compute the value for d2

    double normal(double x);//approximation for normal N(x)

    //function to compute the parameters with given data
    double compute_call_price();
    double compute_put_price();
    double compute_call_delta();
    double compute_put_delta();
    double compute_gamma();
    double compute_vega();
    double compute_call_theta();
    double compute_put_theta();

public:
    //constructor
    bs();
    //rate and div have default argument to be 0
    bs(double s, double k, double vol, double t, double rate, double div);

    //getter functions to get the parameters (price, greeks) of the option
    double call_price();
    double put_price();
    double call_delta();
    double put_delta();
    double gamma();
    double vega();
    double call_theta();
    double put_theta();

};
#endif

