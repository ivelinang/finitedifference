/*
 * Function of the Boundary Condition in a heat PDE u(x, tau) for x = x_right.
 *
 * To be used a functor.
 */

#ifndef GRIGHT_HPP
#define GRIGHT_HPP

#include<math.h>
using namespace std;

class Gright{

    public:
        virtual double operator()(double tau) = 0;
        virtual Gright* clone() const = 0;
};

class Hw8right: public Gright{

    public:
        virtual double operator()(double tau){
            return exp(tau+2.0);
        }

        virtual Gright* clone() const{
            return new Hw8right(*this);
        }
};



#endif
