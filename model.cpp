#include "PrimitiveVar.h"
#include "exactsolver.h"
#include "model.h"
#include "FluxVar.h"
#include "VDC.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <random>

double Model::CFL = 0.8;
double Model::default_NumOfPnt = 251;
double Model::default_NumOfStep = 50;
double Model::default_L = 1;
double Model::default_x_00 = -0.5;

Model::Model(){set(-1,default_NumOfPnt,default_NumOfStep,default_L);}
Model::Model(int choice){set(choice,default_NumOfPnt,default_NumOfStep,default_L);}
Model::Model(int choice, int NumOfStep){set(choice,default_NumOfPnt,NumOfStep,default_L);}
Model::Model(int choice, int _NumOfPnt, int _NumOfStep, double L){set(choice, _NumOfPnt,_NumOfStep,  L);}

double Model::simulate(int solver, bool saveAll){
	// std::string fn = "Data/" + std::to_string(solver)+"S_"+filename + ".txt";
	std::string fn = "Data/p.txt";
    std::ofstream fout(fn.c_str());
    fout << NumOfStep << '\t' << NumOfPnt << std::endl;
    if(saveAll) save_iteration(fout, 0);
if(solver==1){// 1 : GODUNOV, EULER
	LAG = false;
	for(int i = 1; i<=NumOfStep;i++){
		godunov_step(T_act);
		if(saveAll) save_iteration(fout, i);

	}}
if(solver==2){// 2 : GODUNOV, LAGRANGE
	LAG = true;
	for(int i = 1; i<=NumOfStep;i++){
		godunov_lag_step(T_act);
		if(saveAll) save_iteration(fout, i);
	}}
if(solver==3){// 3 : RCM, EULER
	LAG = false;
	double a;
	std::uniform_real_distribution<double> unif(0,1);
   	std::default_random_engine re;
	for(int i = 1; i<=NumOfStep;i++){
		a = unif(re);
		RCM_euler_step(T_act,a);
		if(saveAll) save_iteration(fout, i);
	}}
if(solver==4){// 4 : PSEUDO - RCM, EULER
	LAG = false;
	double* a = vdc_sequence(9,9+NumOfStep+1);
	for(int i = 1; i<=NumOfStep;i++){
		RCM_euler_step(T_act,a[i-1]);
		if(saveAll) save_iteration(fout, T_act);
	}}
if(solver==5){// 5 : RCM, LAGRNAGIAN
	LAG = true;
	double a;
	std::uniform_real_distribution<double> unif(0,1);
   	std::default_random_engine re;
	for(int i = 1; i<=NumOfStep;i++){
		a = unif(re);
		RCM_lag_step(T_act,a);
		if(saveAll) save_iteration(fout, i);
	}}
if(solver==6){// 6 : PSEUDO - RCM, LAGRANGIAN
	LAG = true;
	double* a = vdc_sequence(1,NumOfStep+1);
	for(int i = 0; i<NumOfStep;i++){
		RCM_lag_step(T_act,a[i-1]);
		if(saveAll) save_iteration(fout, i);
	}}
if(solver==7){// 7 : 2-STEPS PSEUDO-RCM, EULER
	LAG = false;
	double* a = vdc_sequence(9,9+2*NumOfStep+2);
	std::vector<PrimitiveVar *>W_inter(NumOfPnt);
	for(int i=0;i<NumOfPnt;i++)
        W_inter[i] = new PrimitiveVar();
	for(int i = 1; i<=NumOfStep;i++){
		RCM_euler_2step(T_act,a[i*2-2], a[2*i-1],W_inter);
		if(saveAll) save_iteration(fout, i);
	}}
return T_act;
}

void Model::godunov_step(double &T){
	// FLUX 
	for(int i=0;i<NumOfPnt;i++)
		ExactSolver::setFlux(*W[i], *W[i+1], *F[i]);
	// PAS DE TEMPS
	double dt = _dt();
	// MISE A JOUR
	for(int i=1; i<NumOfPnt;i++)
		ExactSolver::update(*W[i],*F[i-1]-*F[i], dt/(x[i]-x[i-1]));
	// MISE A JOUR DE LA SOLUTION EXACTE
	T+=dt;
	ExactSolver::exact_sol(W_E,x,T,Wl,Wsl,Wsr,Wr,NumOfPnt);
	W[0]->set(*W_E[0]);
	W[NumOfPnt]->set(*W_E[NumOfPnt]);
}
void Model::godunov_lag_step(double &T){
	// FLUX
	for(int i=0;i<NumOfPnt;i++)
		ExactSolver::setLagFlux(*W[i], *W[i+1], *F[i]);
	// PAS DE TEMPS
	double dt = _dt();
	// MISE A JOUR
	for(int i=1; i<NumOfPnt;i++)
		ExactSolver::lagUpdate(*W[i],*F[i-1]-*F[i], dt/(m[i]-m[i-1]));
	T+=dt;
	// MISE A JOUR DE LA SOLUTION EXACTE
	update_euler_coords(x,m,W, x_00+T*Wl.u);
	ExactSolver::exact_sol(W_E,x,T,Wl,Wsl,Wsr,Wr,NumOfPnt);
}
void Model::RCM_euler_step(double &T, double _a){
	// PAS DE TEMPS
	double dt = _dt();
	PrimitiveVar _Wsl = PrimitiveVar();
	PrimitiveVar _Wsr = PrimitiveVar();
	double p_s,u_s,rho_sL,rho_sR;
	// MISE A JOUR
	PrimitiveVar W_temp = PrimitiveVar();
	if(_a<0.5){
		for(int i=NumOfPnt-2;i>=0;i--){
			p_s = ExactSolver::p_star(*W[i],*W[i+1]);
			u_s = ExactSolver::u_star(p_s,*W[i],*W[i+1]);
		    rho_sL = ExactSolver::rho_star(p_s, *W[i]);
		    rho_sR = ExactSolver::rho_star(p_s, *W[i+1]);
		    _Wsl.set(rho_sL,u_s,p_s);
		    _Wsr.set(rho_sR,u_s,p_s);
		    ExactSolver::Sol_Sample(*W[i],_Wsl,_Wsr,*W[i+1],_a*(x[i+1]-x[i])/dt,W_temp);
		    W[i+1]->set(W_temp);
		}
	}
	else if(_a>=0.5){
		for(int i=1;i<NumOfPnt;i++){
			p_s = ExactSolver::p_star(*W[i],*W[i+1]);
			u_s = ExactSolver::u_star(p_s,*W[i],*W[i+1]);
		    rho_sL = ExactSolver::rho_star(p_s, *W[i]);
		    rho_sR = ExactSolver::rho_star(p_s, *W[i+1]);
		    _Wsl.set(rho_sL,u_s,p_s);
		    _Wsr.set(rho_sR,u_s,p_s);		
		    ExactSolver::Sol_Sample(*W[i],_Wsl,_Wsr,*W[i+1],(_a-1)*(x[i]-x[i-1])/dt,W_temp);
		    W[i]->set(W_temp);
		}
	}
	T+=dt;
	// MISE A JOUR DE LA SOLUTION EXACTE
	ExactSolver::exact_sol(W_E,x,T,Wl,Wsl,Wsr,Wr,NumOfPnt);
}
void Model::RCM_euler_2step(double &T,double a1,double a2,std::vector<PrimitiveVar *> W_inter){
	// PAS DE TEMPS
	double dt = _dt();
	PrimitiveVar _Wsl = PrimitiveVar();
	PrimitiveVar _Wsr = PrimitiveVar();
	double p_s,u_s,rho_sL,rho_sR;
	// MISE A JOUR
	double xL,xR;
	for(int i=0;i<NumOfPnt;i++){
		p_s = ExactSolver::p_star(*W[i],*W[i+1]);
		u_s = ExactSolver::u_star(p_s,*W[i],*W[i+1]);
	    rho_sL = ExactSolver::rho_star(p_s, *W[i]);
	    rho_sR = ExactSolver::rho_star(p_s, *W[i+1]);
	    _Wsl.set(rho_sL,u_s,p_s);
	    _Wsr.set(rho_sR,u_s,p_s);
	    if(i==0) xL = 1.5*x[0]-0.5*x[1];
	    else xL = 0.5*(x[i]+x[i-1]);

	    if(i==NumOfPnt-1) xR = 1.5*x[NumOfPnt-1]-0.5*x[NumOfPnt-2];
	    else xR = 0.5*(x[i]+x[i+1]);
	    ExactSolver::Sol_Sample(*W[i],_Wsl,_Wsr,*W[i+1],((1-a1)*xL+a1*xR-x[i])/dt,*W_inter[i]);
	}
	T+=dt;
	double dX_min=1;
	double S_max = -1;
	for(int i = 1; i<NumOfPnt-1;i++){
		dX_min = std::min(0.5*(x[i+1]-x[i-1]), dX_min);
		S_max = std::max(std::abs(W_inter[i]->u+W_inter[i]->a), S_max);
		S_max = std::max(std::abs(W_inter[i]->u-W_inter[i]->a), S_max);
	}
	S_max = std::max(W_inter[0]->a, S_max);
	dt=CFL*dX_min/S_max;

	for(int i=0;i<NumOfPnt-1;i++){
		p_s = ExactSolver::p_star(*W_inter[i],*W_inter[i+1]);
		u_s = ExactSolver::u_star(p_s,*W_inter[i],*W_inter[i+1]);
	    rho_sL = ExactSolver::rho_star(p_s, *W_inter[i]);
	    rho_sR = ExactSolver::rho_star(p_s, *W_inter[i+1]);
	    _Wsl.set(rho_sL,u_s,p_s);
	    _Wsr.set(rho_sR,u_s,p_s);
	    ExactSolver::Sol_Sample(*W_inter[i],_Wsl,_Wsr,*W_inter[i+1],(a2-0.5)*(x[i+1]-x[i])/dt,*W[i+1]);
	}
	T+=dt;
	// MISE A JOUR DE LA SOLUTION EXACTE
	ExactSolver::exact_sol(W_E,x,T,Wl,Wsl,Wsr,Wr,NumOfPnt);
}
void Model::RCM_lag_step(double &T, double _a){
	// PAS DE TEMPS
	double dt = _dt();
	PrimitiveVar _Wsl = PrimitiveVar();
	PrimitiveVar _Wsr = PrimitiveVar();
	double p_s,u_s,rho_sL,rho_sR;
	// MISE A JOUR
	PrimitiveVar W_temp = PrimitiveVar();
	PrimitiveVar W_temp2 = PrimitiveVar();

	for(int i=1;i<NumOfPnt;i++){
		p_s = ExactSolver::p_star(*W[i-1],*W[i+1]);
		u_s = ExactSolver::u_star(p_s,*W[i-1],*W[i+1]);
	    rho_sL = ExactSolver::rho_star(p_s, *W[i-1]);
	    rho_sR = ExactSolver::rho_star(p_s, *W[i+1]);
	    _Wsl.set(rho_sL,u_s,p_s);
	    _Wsr.set(rho_sR,u_s,p_s);		
	    ExactSolver::Lag_Sol_Sample(*W[i-1],_Wsl,_Wsr,*W[i+1],(_a-0.5)*(m[i]-m[i-1])/dt,W_temp2);
	    if(i>1) W[i-1]->set(W_temp);
	    W_temp.set(W_temp2);
	}
	W[NumOfPnt-1]->set(W_temp);
	
	T+=dt;
	// MISE A JOUR DE LA SOLUTION EXACTE ET DES COORDONNES EULERIENNES
	update_euler_coords(x,m,W, x_00+T*Wl.u);
	ExactSolver::exact_sol(W_E,x,T,Wl,Wsl,Wsr,Wr,NumOfPnt);
}
void Model::save_iteration(std::ofstream &f,int n){
	f << n << std::endl;
	// COORDONNEES
    f<<1.5*x[0]-0.5*x[1]<<std::endl;
    for (int i = 0; i < NumOfPnt-1; i++)
    	f << 0.5*(x[i]+x[i+1]) << std::endl;
    f<<1.5*x[NumOfPnt-1]-0.5*x[NumOfPnt-2]<<std::endl;;
	// VARIABLES
    for (int i = 0; i <= NumOfPnt; i++)
    {
        f << W[i]->rho << '\t' << W[i]->u << '\t' << W[i]->p << '\t' << W[i]->e;
        f << '\t' << W_E[i]->rho << '\t' << W_E[i]->u << '\t' << W_E[i]->p << '\t' << W_E[i]->e << std::endl;
    }
}

void Model::set(int choice, int _NumOfPnt, int _NumOfStep, double L){
// Choix de la condition initiale
	if(choice==1){
		Wl.set(1.0,0.75,1.0);
		Wr.set(0.125,0.0,0.1);
		x_00 = -0.3;
		Tf = 0.2;
	}else if(choice==2){
		Wl.set(1.0,-2.0,0.4);
		Wr.set(1.0,2.0,0.4);
		x_00 = -0.5;
		Tf = 0.15;
	}else if(choice==3){
		Wl.set(1.0,0.0,1000.0);
		Wr.set(1.0,0.0,0.01);
		x_00 = -0.4;
		Tf = 0.012;
	}else if(choice==4){
		Wl.set(5.99924,19.5975,460.894);
		Wr.set(5.99242,-6.19633,46.0950);
		x_00 = -0.4;
		Tf = 0.035;
	}else if(choice==5){
		Wl.set(1.0,-19.59745,1000.0);
		Wr.set(1.0,-19.59745,0.01);
		x_00 = -0.8;
		Tf = 0.012;
	}
	else{
		std::cout<<"Rentrer successivement rhoL, uL, pL, rhoR, uR et pR"<<std::endl;
		std::cin >> Wl.tau >> Wl.u >> Wl.p;
		std::cin >> Wr.tau >> Wr.u >> Wr.p;
		Wl.set(Wl.tau,Wl.u,Wl.p);
		Wr.set(Wr.tau,Wr.u,Wr.p);
		std::cout<<"Rentrer x0"<<std::endl;
		std::cin >> x_00;
		std::cout<<"Rentrer Tf"<<std::endl;
		std::cin >> Tf;
	}
    double ps = ExactSolver::p_star(Wl,Wr);
    double us = ExactSolver::u_star(ps,Wl,Wr);
    double rhoL = ExactSolver::rho_star(ps,Wl);
    double rhoR = ExactSolver::rho_star(ps,Wr);
    Wsl.set(rhoL,us,ps);
    Wsr.set(rhoR,us,ps);
// Initialisation des coordonnées et des conditions initiales
    T_act=0;
    NumOfPnt = _NumOfPnt;
    NumOfStep = _NumOfStep;
    W = std::vector<PrimitiveVar *>(NumOfPnt+1);
    W_E = std::vector<PrimitiveVar *>(NumOfPnt+1);
    F =std::vector<FluxVar *>(NumOfPnt);

	for(int i=0;i<=NumOfPnt;i++)
		W_E[i] = new PrimitiveVar();
	
	for(int i=0; i<NumOfPnt;i++)
		F[i] = new FluxVar();


	double dX = L/(NumOfPnt-1);
    m = std::vector<double>(NumOfPnt);
	x = std::vector<double>(NumOfPnt);
	x[0] = x_00;
	m[0] = -x_00*Wl.rho;
	W[0] = new PrimitiveVar(Wl);W_E[0] = new PrimitiveVar(Wl);
	for(int i=1;i<NumOfPnt;i++)
	{
		x[i]=x[i-1]+dX;
    	if(0.5*(x[i]+x[i-1])<0){W[i] = new PrimitiveVar(Wl);W_E[i] = new PrimitiveVar(Wl);}
    	else{W[i] = new PrimitiveVar(Wr);W_E[i] = new PrimitiveVar(Wr);}
    	m[i]=m[i-1]+W[i]->rho*dX;
	}
	W[NumOfPnt] = new PrimitiveVar(Wr);W_E[NumOfPnt] = new PrimitiveVar(Wr);


	//x_00 = _x_00;
// Nom pour le fichier
filename = std::to_string(choice)+"_"+  std::to_string(_NumOfPnt)+"_"+  std::to_string(_NumOfStep)+"_"+  std::to_string(L)+"_"+ std::to_string(x_00);
}
double Model::_dt(){
	double dxmin, Smax;
	if(LAG){
		dxmin = 10000000;
		Smax = -1;
		for(int i = 1; i<NumOfPnt;i++){
			dxmin = std::min(dxmin,m[i]-m[i-1]);
			Smax = std::max(Smax, W[i]->a_lag);
		}
	}
	else{
		dxmin = 1;
		Smax = -1;
		for(int i = 1; i<NumOfPnt;i++){
			dxmin = std::min(dxmin,x[i]-x[i-1]);
			Smax = std::max(Smax, std::abs(W[i]->u+W[i]->a));
			Smax = std::max(Smax, std::abs(W[i]->u-W[i]->a));
		}
	}
	return CFL*dxmin/Smax;
}

void Model::update_euler_coords(std::vector<double> &x,const std::vector<double> &m, const std::vector<PrimitiveVar *> w, double x_0){
    x[0]=x_0;
    for(int i = 1; i<NumOfPnt;i++)
        x[i]=x[i-1]+(m[i]-m[i-1])*w[i]->tau;
}

void Model::update_lag_coords(const std::vector<double> &x,std::vector<double> &m, const std::vector<PrimitiveVar *> w){
    double m_tot =0;
    for(int i=1;i<NumOfPnt;i++)
    	m_tot+=W[i]->tau*(x[i]-x[i-1]);
    m[0]=-m_tot/2;
    for(int i=1; i<NumOfPnt;i++)
    	m[i] = m[i-1]+W[i]->tau*(x[i]-x[i-1]);
}
