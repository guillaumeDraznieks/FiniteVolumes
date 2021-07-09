#include "PrimitiveVar.h"
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



double PrimitiveVar::TOL = 1e-6;
double PrimitiveVar::G0 = 1.4;
double PrimitiveVar::G1 = 0.5*(G0 - 1) / G0;
double PrimitiveVar::G2 = 0.5*(G0 + 1) / G0;
double PrimitiveVar::G3 = 2 * G0 / (G0 - 1);
double PrimitiveVar::G4 = 2 / (G0 - 1);
double PrimitiveVar::G5 = 2 / (G0 + 1);
double PrimitiveVar::G6 = (G0 - 1) / (G0 + 1);
double PrimitiveVar::G7 = (G0 - 1) / 2;
double PrimitiveVar::G8 = G0 - 1;
double PrimitiveVar::G9 = -G2;
double PrimitiveVar::G10 = -0.5*(G0 + 1) / pow(G0, 2);
double PrimitiveVar::G11 = -0.5*(3 * G0 + 1) / G0;
double PrimitiveVar::G12 = 1 / G0;


PrimitiveVar::PrimitiveVar(double density, double velocity, double pressure)
{
    set(density, velocity, pressure);
}
PrimitiveVar::PrimitiveVar(){
    set(0,0,0);
}
PrimitiveVar::PrimitiveVar(const PrimitiveVar &W){
    set(W);
}


void PrimitiveVar::set(double density, double velocity, double pressure)
{
    rho = density;
    u = velocity;
    p = pressure;
    tau = 1/density;
    a = sound_speed(p, rho);
    e = internal_energy(p, rho);
    a_lag = lag_sound_speed(p,tau);
}
void PrimitiveVar::lag_set(double _tau, double velocity, double pressure)
{
    rho = 1/_tau;
    u = velocity;
    p = pressure;
    tau = _tau;
    a = sound_speed(p, rho);
    a_lag = lag_sound_speed(p,tau);
    e = internal_energy(p, rho);
}
void PrimitiveVar::set(const PrimitiveVar &W){
    rho = W.rho;
    u=W.u;
    p=W.p;
    a=W.a;
    e=W.e;
    tau=W.tau;
    a_lag = W.a_lag;
}




inline double PrimitiveVar::E(double density, double velocity, double pressure)
{
    return density * (internal_energy(pressure, density) + kinetic_energy(velocity));
}
inline double PrimitiveVar::E(const PrimitiveVar &W)
{
    return W.rho * (W.e + kinetic_energy(W.u));
}
inline double PrimitiveVar::kinetic_energy(double u)
{
    return 0.5 * pow(u, 2);
}
inline double PrimitiveVar::sound_speed(double p, double rho)
{
    return sqrt(G0 * p / rho);
}
inline double PrimitiveVar::lag_sound_speed(double p, double tau)
{
    return sqrt(G0 * p / tau);
}
inline double PrimitiveVar::internal_energy(double p, double rho)
{
    return p / (G8 * rho);
}






