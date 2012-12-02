/*
 * Function of the Boundary Condition in a heat PDE u(x, tau) for x = x_left.
 *
 * To be used a functor.
 */

#ifndef GLEFT_HPP
#define GLEFT_HPP

#include<math.h>
using namespace std;

class Gleft{

    public:
        virtual double operator()(double tau) = 0;
        virtual Gleft* clone() const = 0;
};

class Hw8left: public Gleft{

    public:
        virtual double operator()(double tau){
            return exp(tau-2.0);
        }

        virtual Gleft* clone() const{
            return new Hw8left(*this);
        }
};



#endif
