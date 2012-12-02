/*
 * Definition of the heat pde solver interface
 *
 * Yongyi Ye
 */

#include"heatpdesolver.hpp"


HeatPdeSolver::HeatPdeSolver(double xleft_, double xright_, double taufinal_,
                               const Gleft &gleft_, const Gright &gright_, const Ftau &f_):
                            xleft(xleft_), xright(xright_), taufinal(taufinal_),
                              gleft(gleft_), gright(gright_), f(f_){}


HeatPdeSolver::HeatPdeSolver(const HeatPdeSolver &input):xleft(input.xleft), xright(input.xright), taufinal(input.taufinal),
                                                         gleft(input.gleft), gright(input.gright), f(input.f){}

HeatPdeSolver::~HeatPdeSolver(){}

HeatPdeSolver& HeatPdeSolver::operator= (const HeatPdeSolver &input){
    xleft = input.xleft;
    xright = input.xright;
    taufinal = input.taufinal;
    gleft = input.gleft;
    gright = input.gright;
    f = input.f;
}
