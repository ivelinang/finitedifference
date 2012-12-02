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



#endif
