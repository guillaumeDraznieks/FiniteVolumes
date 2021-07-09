#include "exactsolver.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

using namespace std;

double ExactSolver::TOL = 1e-6;
double ExactSolver::G0 = 1.4;
double ExactSolver::G1 = 0.5*(G0 - 1) / G0;
double ExactSolver::G2 = 0.5*(G0 + 1) / G0;
double ExactSolver::G3 = 2 * G0 / (G0 - 1);
double ExactSolver::G4 = 2 / (G0 - 1);
double ExactSolver::G5 = 2 / (G0 + 1);
double ExactSolver::G6 = (G0 - 1) / (G0 + 1);
double ExactSolver::G7 = (G0 - 1) / 2;
double ExactSolver::G8 = G0 - 1;
double ExactSolver::G9 = -G2;
double ExactSolver::G10 = -0.5*(G0 + 1) / pow(G0, 2);
double ExactSolver::G11 = -0.5*(3 * G0 + 1) / G0;
double ExactSolver::G12 = 1 / G0;
double ExactSolver::G13 = 2/(G0+1);



inline double kinetic_energy(double u)
{
    return 0.5 * pow(u, 2);
}



void ExactSolver::setFlux(const PrimitiveVar &Wl, const PrimitiveVar &Wr,FluxVar &F){
	double p_s = p_star(Wl,Wr);
	double u_s = u_star(p_s,Wl,Wr);
    double rho_sL = rho_star(p_s, Wl);
    double rho_sR = rho_star(p_s, Wr);

    //finding the flux

    if (0 < u_s)
        {
            //At the left of the contact discontinuity
            if (p_s > Wl.p)
            {
                //Left shock
                double S_L = Wl.u - Wl.a * sqrt(G2 * p_s / Wl.p + G1);
                if (0 < S_L)
                	F.set(Wl.rho, Wl.u, Wl.p);
                else
                	F.set(rho_sL, u_s, p_s);
            }
            else
            {
                //Left fan
                double S_HL = Wl.u - Wl.a;
                if (0 < S_HL)
                    F.set(Wl.rho, Wl.u, Wl.p);
                else
                {
                    double S_TL = u_s - sound_speed(p_s, rho_sL);
                    if (0 > S_TL)
                        F.set(rho_sL, u_s, p_s);
                    else{
						double density = Wl.rho * pow(G5 + G6 / Wl.a * Wl.u, G4);
						double velocity = G5 * (Wl.a + G7 * Wl.u);
						double pressure = Wl.p * pow(G5 + G6 / Wl.a * Wl.u , G3);
						F.set(density, velocity, pressure);
                    }
                }
            }
        }
        else
        {
            //At the right of the contact discontinuity
            if (p_s > Wr.p)
            {
                //Right shock
                double S_R = Wr.u + Wr.a * sqrt(G2 * p_s / Wr.p + G1);
                if (0 < S_R)
                    F.set(rho_sR, u_s, p_s);
                else
                    F.set(Wr.rho, Wr.u, Wr.p);
            }
            else
            {
                //Right fan
                double S_HR = Wr.u + Wr.a;
                if (0 > S_HR)
                    F.set(Wr.rho, Wr.u, Wr.p);
                else
                {
                    double S_TR = u_s + sound_speed(p_s, rho_sR);
                    if (0 < S_TR)
                        F.set(rho_sR, u_s, p_s);
                    else{
                        double density = Wr.rho * pow(G5 - G6 / Wr.a * Wr.u, G4);
					    double velocity = G5 * (-Wr.a + G7 * Wr.u);
					    double pressure = Wr.p * pow(G5 - G6 / Wr.a * Wr.u, G3);
					    F.set(density, velocity, pressure);
                    }
                }
            }
        }
}
void ExactSolver::setLagFlux(const PrimitiveVar &Wl, const PrimitiveVar &Wr,FluxVar &F){
	double p_s = p_star(Wl,Wr);
	double u_s = u_star(p_s,Wl,Wr);
	F.lag_set(u_s,p_s);
}
void ExactSolver::update(PrimitiveVar &W,const FluxVar &dF,double dt){
    // variables conservatives
    double _rho = W.rho;
    double _rho_u = _rho*W.u;
    double _E = _rho*(kinetic_energy(W.u)+W.e);

    //on update
    _rho += dt*dF.f[0];
    _rho_u += dt*dF.f[1];
    _E += dt*dF.f[2];

    //on update les variables primaires
    double _u = _rho_u/_rho;
    W.set(_rho, _u, 0.4*(_E-_rho*kinetic_energy(_u)));
}
void ExactSolver::lagUpdate(PrimitiveVar &W,const FluxVar &dF,double dt){
    // variables conservatives
    double _tau = W.tau;
    double _u = W.u;
    double _E = kinetic_energy(W.u)+W.e;

    //on update
    _tau += dt*dF.f[0];
    _u += dt*dF.f[1];
    _E += dt*dF.f[2];
    
    W.lag_set(_tau, _u, 0.4*(_E-kinetic_energy(_u))/_tau);
}
void ExactSolver::exact_sol(std::vector<PrimitiveVar *> &w, const std::vector<double> &x, double t, const PrimitiveVar &Wl, const PrimitiveVar &Wsl, const PrimitiveVar &Wsr, const PrimitiveVar &Wr, int NumOfPnt){
	ExactSolver::Sol_Sample(Wl, Wsl, Wsr, Wr, x[0]/t,*w[0]);
	ExactSolver::Sol_Sample(Wl, Wsl, Wsr, Wr, x[NumOfPnt-1]/t,*w[NumOfPnt]);
	for(int i = 1; i<NumOfPnt;i++)
    	ExactSolver::Sol_Sample(Wl, Wsl, Wsr, Wr, 0.5*(x[i]+x[i-1])/t,*w[i]);
}
void ExactSolver::Sol_Sample(const PrimitiveVar &Wl, const PrimitiveVar &Wsl, const PrimitiveVar &Wsr, const PrimitiveVar &Wr, double S, PrimitiveVar &ans)
{
	double u_star = Wsl.u;
	double S_L = Wl.u - Wl.a * sqrt(G2 * Wsl.p / Wl.p + G1);
	double S_HL = Wl.u - Wl.a;
	double S_TL = Wsl.u - Wsl.a;
	double S_R = Wr.u + Wr.a * sqrt(G2 * Wsr.p / Wr.p + G1);
	double S_HR = Wr.u + Wr.a;
	double S_TR = Wsr.u + Wsr.a;
	

	if (Wsl.p > Wl.p)
	{
		//Left Shock
		if (Wsr.p > Wr.p)
		{
			//Right Shock
			if (S < S_L)
				ans.set(Wl);
			else if (S < u_star)
				ans.set(Wsl);
			else if (S < S_R)
				ans.set(Wsr);
			else
				ans.set(Wr);
		}
		else
		{
			//Right Rarefaction
			if (S < S_L)
				ans.set(Wl);
			else if (S < u_star)
				ans.set(Wsl);
			else if (S < S_TR)
				ans.set(Wsr);
			else if (S < S_HR)
				Calc_WfanR(Wr, S, ans);
			else
				ans.set(Wr);
		}
	}
	else
	{
		//Left Rarefaction
		if (Wsr.p > Wr.p)
		{
			//Right Shock
			if (S < S_HL)
				ans.set(Wl);
			else if (S < S_TL)
				Calc_WfanL(Wl, S, ans);
			else if (S < u_star)
				ans.set(Wsl);
			else if (S < S_R)
				ans.set(Wsr);
			else
				ans.set(Wr);
		}
		else
		{
			//Right Rarefaction
			if (S < S_HL)
				ans.set(Wl);
			else if (S < S_TL)
				Calc_WfanL(Wl, S, ans);
			else if (S < u_star)
				ans.set(Wsl);
			else if (S < S_TR)
				ans.set(Wsr);
			else if (S < S_HR)
				Calc_WfanR(Wr, S, ans);
			else
				ans.set(Wr);
		}
	}
}
void ExactSolver::Lag_Sol_Sample(const PrimitiveVar &Wl, const PrimitiveVar &Wsl, const PrimitiveVar &Wsr, const PrimitiveVar &Wr, double S, PrimitiveVar &ans)
{
	double S_L = (Wsl.p-Wl.p)/(Wsl.u-Wl.u);
	double S_HL = - Wl.a_lag;
	double S_TL = -	Wsl.a_lag;
	double S_R = (Wsr.p-Wr.p)/(Wsr.u-Wr.u);
	double S_HR = Wr.a_lag;
	double S_TR = Wsr.a_lag;
	

	if (Wsl.p > Wl.p)
	{
		//Left Shock
		if (Wsr.p > Wr.p)
		{
			//Right Shock
			if (S < S_L)
				ans.set(Wl);
			else if (S < 0)
				ans.set(Wsl);
			else if (S < S_R)
				ans.set(Wsr);
			else
				ans.set(Wr);
		}
		else
		{
			//Right Rarefaction
			if (S < S_L)
				ans.set(Wl);
			else if (S < 0)
				ans.set(Wsl);
			else if (S < S_TR)
				ans.set(Wsr);
			else if (S < S_HR)
				Lag_Calc_WfanR(Wr, S, ans);
			else
				ans.set(Wr);
		}
	}
	else
	{
		//Left Rarefaction
		if (Wsr.p > Wr.p)
		{
			//Right Shock
			if (S < S_HL)
				ans.set(Wl);
			else if (S < S_TL)
				Lag_Calc_WfanL(Wl, S, ans);
			else if (S < 0)
				ans.set(Wsl);
			else if (S < S_R)
				ans.set(Wsr);
			else
				ans.set(Wr);
		}
		else
		{
			//Right Rarefaction
			if (S < S_HL)
				ans.set(Wl);
			else if (S < S_TL)
				Lag_Calc_WfanL(Wl, S, ans);
			else if (S < 0)
				ans.set(Wsl);
			else if (S < S_TR)
				ans.set(Wsr);
			else if (S < S_HR)
				Lag_Calc_WfanR(Wr, S, ans);
			else
				ans.set(Wr);
		}
	}
}
void ExactSolver::Calc_WfanL(const PrimitiveVar &Wk, double S, PrimitiveVar &ans)
{
	double rho = Wk.rho * pow(G5 + G6 / Wk.a *(Wk.u - S), G4);
	double u = G5 * (Wk.a + G7 * Wk.u + S);
	double p = Wk.p * pow(G5 + G6 / Wk.a * (Wk.u - S), G3);
	ans.set(rho,u,p);
}
void ExactSolver::Calc_WfanR(const PrimitiveVar &Wk, double S, PrimitiveVar &ans)
{
	double rho = Wk.rho * pow(G5 - G6 / Wk.a *(Wk.u - S), G4);
	double u = G5 * (-Wk.a + G7 * Wk.u + S);
	double p = Wk.p * pow(G5 - G6 / Wk.a * (Wk.u - S), G3);
	ans.set(rho,u,p);
}
void ExactSolver::Lag_Calc_WfanL(const PrimitiveVar &Wk, double S, PrimitiveVar &ans)
{
	double i = S/Wk.a_lag;
	double rho = Wk.rho * pow(abs(i), G5);
	double u = Wk.u + G4*Wk.a*(1+i*Wk.rho/rho);
	double p = Wk.p * pow(abs(i),G13);
	ans.set(rho,u,p);
}
void ExactSolver::Lag_Calc_WfanR(const PrimitiveVar &Wk, double S, PrimitiveVar &ans)
{
	double i = S/Wk.a_lag;
	double rho = Wk.rho * pow(abs(i), G5);
	double u = Wk.u + G4*Wk.a*(i*Wk.rho/rho-1);
	double p = Wk.p * pow(abs(i),G13);
	ans.set(rho,u,p);
}
inline double ExactSolver::rho_star(double p, const PrimitiveVar &W)
{
    const double t = p / W.p;

    if (t > 1.0)
        return W.rho * ((t + G6) / (G6 * t + 1));
    else
        return W.rho * pow(t, G12);
}
inline double ExactSolver::p_star(const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
    const double du = Wr.u - Wl.u;

    //Select initial pressure
    const double f_min = f(min(Wl.p, Wr.p), Wl, Wr);
    const double f_max = f(max(Wl.p, Wr.p), Wl, Wr);

    double p_m = (Wl.p + Wr.p) / 2;
    p_m = max(TOL, p_m);

    double p_pv = p_m - du * (Wl.rho + Wr.rho) * (Wl.a + Wr.a) / 8;
    p_pv = max(TOL, p_pv);

    const double gL = gk(p_pv, Wl);
    const double gR = gk(p_pv, Wr);
    double p_ts = (gL * Wl.p + gR * Wr.p - du) / (gL + gR);
    p_ts = max(TOL, p_ts);

    double p_tr = pow((Wl.a + Wr.a - G7 * du) / (Wl.a / pow(Wl.p, G1) + Wr.a / pow(Wr.p, G1)), G3);
    p_tr = max(TOL, p_tr);

    double p0 = p_m;
    if (f_min < 0 && f_max < 0)
        p0 = p_ts;
    else if (f_min > 0 && f_max > 0)
        p0 = p_tr;
    else
        p0 = p_pv;

    //Solve
    int iter_cnt = 0;
    double CHA = 1.0;
    while (CHA > TOL)
    {
        ++iter_cnt;

        double fder = df(p0, Wl, Wr);
        if (fder == 0)
            throw "Zero derivative!";

        double fval = f(p0, Wl, Wr);
        double fder2 = ddf(p0, Wl, Wr);
        double p = p0 - fval * fder / (pow(fder, 2) - 0.5 * fval * fder2);
        if (p < 0)
        {
            p0 = TOL;
            break;
        }

        CHA = abs(2 * (p - p0) / (p + p0));
        p0 = p;
    }

    return p0;
}
inline double ExactSolver::u_star(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
    return (Wl.u + Wr.u) / 2 + (fk(p, Wr) - fk(p, Wl)) / 2;
}
inline double ExactSolver::sound_speed(double p, double rho)
{
	return sqrt(G0*p / rho);
}
inline double ExactSolver::Ak(const PrimitiveVar &Wk)
{
	return G5 / Wk.rho;
}
inline double ExactSolver::Bk(const PrimitiveVar &Wk)
{
	return G6 * Wk.p;
}
inline double ExactSolver::gk(double p, const PrimitiveVar &Wk)
{
	return sqrt(Ak(Wk) / (p + Bk(Wk)));
}
inline double ExactSolver::fk(double p, const PrimitiveVar &Wk)
{
	if (p > Wk.p)
	{
		//shock
		return (p - Wk.p)*gk(p, Wk);
	}
	else
	{
		//rarefraction
		return G4 * Wk.a*(pow(p / Wk.p, G1) - 1);
	}
}
inline double ExactSolver::f(double p, const PrimitiveVar &Wl,const PrimitiveVar &Wr)
{
	return fk(p, Wl) + fk(p, Wr) + (Wr.u - Wl.u);
}

inline double ExactSolver::dfk(double p,const PrimitiveVar &Wk)
{
	if (p > Wk.p)
	{
		//shock
		double A = Ak(Wk);
		double B = Bk(Wk);
		return sqrt(A / (p + B))*(1 - 0.5*(p - Wk.p) / (B + p));
	}
	else
	{
		//rarefraction
		return pow(p / Wk.p, G9) / (Wk.rho*Wk.a);
	}
}
inline double ExactSolver::df(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
	return dfk(p, Wl) + dfk(p, Wr);
}
inline double ExactSolver::ddfk(double p, const PrimitiveVar &Wk)
{
	if (p > Wk.p)
	{
		//shock
		double A = Ak(Wk);
		double B = Bk(Wk);
		return -0.25 * sqrt(A / (p + B)) * (4 * B + 3 * p + Wk.p) / pow(B + p, 2);
	}
	else
	{
		//rarefraction
		return G10 * Wk.a / pow(Wk.p, 2) * pow(p / Wk.p, G11);
	}
}
inline double ExactSolver::ddf(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
	return ddfk(p, Wl) + ddfk(p, Wr);
}


