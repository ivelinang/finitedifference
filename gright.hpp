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
        virtual ~Gright();
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

class BsEuropeanPutRight: public Gright{
/* the boundary condition gright for the underlying heat pde for solving the black-scholes pde of a EUROPEAN PUT */

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

    public:
        BsEuropeanPutRight(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_right boundary value for the heat pde such that x_right > x */
        virtual double get_x_right();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gright* clone() const;
};

class BsEuropeanCallRight: public Gright{
/* the boundary condition gright for the underlying heat pde for solving the black-scholes pde of a EUROPEAN CALL*/

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
        BsEuropeanCallRight(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_right boundary value for the heat pde such that x_right > x */
        virtual double get_x_right();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gright* clone() const;
};

class BsAmericanPutRight: public Gright{
/* the boundary condition gright for the underlying heat pde for solving the black-scholes pde of a AMERICAN PUT */

    private:
        double s;   // asset price at time 0
        double k;   // options strike
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

    public:
        BsAmericanPutRight(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_right boundary value for the heat pde such that x_right > x */
        virtual double get_x_right();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gright* clone() const;
};

class BsAmericanCallRight: public Gright{
/* the boundary condition gright for the underlying heat pde for solving the black-scholes pde of a AMERICAN CALL*/

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
        BsAmericanCallRight(double s_, double k_, double vol_, double t_, double r_, double q_);

        /* get the x_right boundary value for the heat pde such that x_right > x */
        virtual double get_x_right();

        /* the boundary condition equation */
        virtual double operator()(double tau);

        /* virtual copy constructor */
        virtual Gright* clone() const;
};

#endif
