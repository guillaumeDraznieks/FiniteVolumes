#ifndef MODEL_H
#define MODEL_H

#include "FluxVar.h"
#include "PrimitiveVar.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>

class Model
{
  public:

    static double CFL;
    static double default_NumOfPnt;
    static double default_NumOfStep;
    static double default_L;
    static double default_x_00;
  // QUANTITES UTILES
    int NumOfPnt;
    int NumOfStep;
    double x_00;double T_act;
    double Tf;
    bool LAG;
  // FICHIER
    std::string filename;
  // SOLUTION EXACT
    PrimitiveVar Wl,Wsl,Wsr,Wr;
  // VECTEURS DE TRAVAIL


    Model();
    Model(int choice);
    Model(int choice, int NumOfStep);
    Model(int choice, int NumOfPnt, int NumOfStep, double L);

    double simulate(int solver, bool saveAll);
    // 1 : EXACT GODUNOV EULER
    // 2 : EXACT GODUNOV LAGRANGIAN
    // 3 : RCM, EULER
    // 4 : PSEUDO-RCM, EULER
    // 5 : RCM, LAGRANGIAN
    // 6 : PSEUDO-RCM, LAGRANGIAN
    // 7 : 2-STEPS PSEUDO-RCM, EULER

    double Lp_error(int p);

private:

    std::vector<double> x;
    std::vector<double> m;
    std::vector<PrimitiveVar *> W;
    std::vector<PrimitiveVar *> W_E;
    std::vector<FluxVar *> F;

    void set(int choice, int _NumOfPnt, int _NumOfStep, double L);
    void update_euler_coords(std::vector<double> &x,const std::vector<double> &m, const std::vector<PrimitiveVar *> w, double x_0);
    void update_lag_coords(const std::vector<double> &x,std::vector<double> &m, const std::vector<PrimitiveVar *> w);
    void godunov_step(double &T);
    void godunov_lag_step(double &T); 
    void RCM_euler_step(double &T,double a);//a est un nombre aléatoire entre 0 et 1 à fournir
    void RCM_euler_2step(double &T,double a1,double a2,std::vector<PrimitiveVar *> W_inter);
    void RCM_lag_step(double &T,double a);
    double _dt();
    void save_iteration(std::ofstream &f,int n);

};

#endif //MODEL_H