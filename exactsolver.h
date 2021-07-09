#ifndef EXACTSOLVER_H
#define EXACTSOLVER_H

#include "PrimitiveVar.h"
#include "FluxVar.h"

#include <vector>

class ExactSolver
{
public:


    static double TOL;
    static double G0,G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12,G13;

    static void setFlux(const PrimitiveVar &Wl, const PrimitiveVar &Wr,FluxVar &F);
    static void setLagFlux(const PrimitiveVar &Wl, const PrimitiveVar &Wr,FluxVar &F);
    static void update(PrimitiveVar &W, const FluxVar &dF,double dt);
    static void lagUpdate(PrimitiveVar &W, const FluxVar &dF,double dt);
    static void exact_sol(std::vector<PrimitiveVar *> &w, const std::vector<double> &x, double t, const PrimitiveVar &Wl, const PrimitiveVar &Wsl, const PrimitiveVar &Wsr, const PrimitiveVar &Wr, int NumOfPnt);
    
    static double rho_star(double p, const PrimitiveVar &W);
    static double p_star(const PrimitiveVar &Wl, const PrimitiveVar &Wr);
    static double u_star(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr);
    
    static void Sol_Sample(const PrimitiveVar &Wl,const  PrimitiveVar &Wsl,const  PrimitiveVar &Wsr,const  PrimitiveVar &Wr, double S, PrimitiveVar &ans);
    static void Lag_Sol_Sample(const PrimitiveVar &Wl, const PrimitiveVar &Wsl, const PrimitiveVar &Wsr, const PrimitiveVar &Wr, double S, PrimitiveVar &ans);

private:

    static void Calc_WfanL(const PrimitiveVar &Wk, double S, PrimitiveVar &ans);
    static void Calc_WfanR(const PrimitiveVar &Wk, double S, PrimitiveVar &ans);    
    static void Lag_Calc_WfanL(const PrimitiveVar &Wk, double S, PrimitiveVar &ans);
    static void Lag_Calc_WfanR(const PrimitiveVar &Wk, double S, PrimitiveVar &ans);        

    static double sound_speed(double p, double rho);
    static double Ak(const PrimitiveVar &Wk);
    static double Bk(const PrimitiveVar &Wk);
    static double gk(double p, const PrimitiveVar &Wk);
    static double fk(double p, const PrimitiveVar &Wk);
    static double f(double p, const PrimitiveVar &Wl,const PrimitiveVar &Wr);
    static double dfk(double p,const PrimitiveVar &Wk);
    static double df(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr);
    static double ddfk(double p, const PrimitiveVar &Wk);
    static double ddf(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr);
};

#endif // EXACTSOLVER_H