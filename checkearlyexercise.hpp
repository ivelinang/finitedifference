/*
 *  Declaration:
 *
 *  CheckEarlyExercise is the class that will be used by the tailored version of Forward Euler and
 *  Crank-Nicolson that are changed to be used for the American options with early exercise.
 *  It will be used by the heat pde solvers to check for knots of early exercise when solving the
 *  heat pde of x and tau.
 *
 *  Yongyi Ye
 */

#ifndef CHECKEARLYEXERCISE_HPP
#define CHECKEARLYEXERCISE_HPP

class CheckEarlyExercise{

    protected:
        double s;   // stock price
        double k;   // options strike price
        double vol; // volatility
        double t;   // maturity
        double r;   // risk-free interest
        double q;   // continuous dividend rate

        /* constants involved for computing the early exercise condition */
        double x_left;
        double x_domain;
        double tau_domain;
        double a;
        double b;

    public:
        CheckEarlyExercise(double s_, double k_, double vol_, double t_, double r_, double q_);
        CheckEarlyExercise(const CheckEarlyExercise &input);
        virtual ~CheckEarlyExercise();
        virtual CheckEarlyExercise& operator= (const CheckEarlyExercise &input);

        virtual double check_early_exercise(int N, int M, int n, int m) = 0;
        /* function to check for early exercise
         * input: N, M - the dimension of the finite difference grid that the computational domain is divided into
         *               when solving the heat pde
         *        n, m - the current position of x and tau in the finite difference grid
         */

        virtual CheckEarlyExercise* clone() const = 0;
        /* virtual copy constructor */
};

class CheckAmericanPut: public CheckEarlyExercise{

    public:
        CheckAmericanPut(double s_, double k_, double vol_, double t_, double r_, double q_);
        CheckAmericanPut(const CheckAmericanPut &input);
        virtual ~CheckAmericanPut();
        virtual CheckAmericanPut& operator= (const CheckAmericanPut &input);

        virtual double check_early_exercise(int N, int M, int n, int m);
        virtual CheckEarlyExercise* clone() const;
};

class CheckAmericanCall: public CheckEarlyExercise{

    public:
        CheckAmericanCall(double s_, double k_, double vol_, double t_, double r_, double q_);
        CheckAmericanCall(const CheckAmericanCall &input);
        virtual ~CheckAmericanCall();
        virtual CheckAmericanCall& operator= (const CheckAmericanCall &input);

        virtual double check_early_exercise(int N, int M, int n, int m);
        virtual CheckEarlyExercise* clone() const;
};

#endif
