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
        virtual ~Gleft();
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

class BsEuropeanPutLeft: public Gleft{
/* the boundary condition gleft for the underlying heat pde for solving the black-scholes pde of a EUROPEAN PUT */

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

        /* constants involved to transform BS PDE to Heat PDE */
        double a;
        double b;

    public:
        BsEuropeanPutLeft(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_left boundary value for the heat pde such that x_left < x */
        virtual double get_x_left();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gleft* clone() const;
};

class BsEuropeanCallLeft: public Gleft{
/* the boundary condition gleft for the underlying heat pde for solving the black-scholes pde of a EUROPEAN CALL*/

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

    public:
        BsEuropeanCallLeft(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_left boundary value for the heat pde such that x_left < x */
        virtual double get_x_left();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gleft* clone() const;
};

class BsAmericanPutLeft: public Gleft{
/* the boundary condition gleft for the underlying heat pde for solving the black-scholes pde of a AMERICAN PUT */

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

        /* constants involved to transform BS PDE to Heat PDE */
        double a;
        double b;

    public:
        BsAmericanPutLeft(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_left boundary value for the heat pde such that x_left < x */
        virtual double get_x_left();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gleft* clone() const;
};

class BsAmericanCallLeft: public Gleft{
/* the boundary condition gleft for the underlying heat pde for solving the black-scholes pde of a AMERICAN CALL*/

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

    public:
        BsAmericanCallLeft(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_left boundary value for the heat pde such that x_left < x */
        virtual double get_x_left();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gleft* clone() const;
};


#endif
