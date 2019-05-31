/*
 * Implementation of the Heat PDE solver with Early Exercise
 *
 * Yongyi Ye
 */

#include"heatpdesolver.hpp"
#include"../blackscholes/checkearlyexercise.hpp"

EarlyExerciseSolver::EarlyExerciseSolver(double xleft_, double xright_, double taufinal_,
                                            const Gleft &gleft_, const Gright &gright_, const Ftau &f_,
                                            const CheckEarlyExercise &checker_):
                                HeatPdeSolver(xleft_, xright_, taufinal_, gleft_, gright_, f_), checker(checker_){}



EarlyExerciseSolver::EarlyExerciseSolver(const EarlyExerciseSolver &input): HeatPdeSolver(input),
                                                                            checker(input.checker){}

EarlyExerciseSolver::~EarlyExerciseSolver(){}

EarlyExerciseSolver& EarlyExerciseSolver::operator= (const EarlyExerciseSolver &input){
    HeatPdeSolver::operator=(input);
    checker = input.checker;

	return *this;
}
