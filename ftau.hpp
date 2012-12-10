/*
 * Function of the Boundary Condition in a heat PDE u(x, tau) for tau = tau_final
 *
 * To be used a functor.
 */

#ifndef FTAU_HPP
#define FTAU_HPP

#include<math.h>
using namespace std;

class Ftau{

    public:
        virtual ~Ftau();
        virtual double operator()(double x) = 0;
        virtual Ftau* clone() const = 0;
};

class Hw8tau: public Ftau{

    public:
        virtual double operator()(double x){
            return exp(x);
        }

        virtual Ftau* clone() const{
            return new Hw8tau(*this);
        }
};

class BsPutTau: public Ftau{
/* the boundary condition f for tau = 0 for the underlying heat pde for solving the black-scholes pde of a PUT
 * same for European and American
 */

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

    public:
        BsPutTau(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the tau_final for the heat pde such that tau < tau_final */
        virtual double get_tau_final();

        /* the boundary condition equation */
        virtual double operator()(double x);

        /* virtual copy constructor */
        virtual Ftau* clone() const;
};

class BsCallTau: public Ftau{
/* the boundary condition f for tau = 0 for the underlying heat pde for solving the black-scholes pde of a CALL
 * same for European and American
 */


    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

    public:
        BsCallTau(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the tau_final for the heat pde such that tau < tau_final */
        virtual double get_tau_final();

        /* the boundary condition equation */
        virtual double operator()(double x);

        /* virtual copy constructor */
        virtual Ftau* clone() const;
};


#endif
