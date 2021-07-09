#include "FluxVar.h"
#include <cmath>
#include <iostream>

using namespace std;


FluxVar::FluxVar(double a, double b, double c){
  f[0]=a;f[1]=b;f[2]=c;
}

FluxVar::FluxVar(){
  f[0]=0;f[1]=0;f[2]=0;
}
void FluxVar::lag_set(double u, double p){
  f[0] = -u;
  f[1] = p;
  f[2] = u * p;
}

void FluxVar::set(double rho, double u, double p){
  f[0] = rho * u;
  f[1] = rho * pow(u, 2) + p;
  double E = rho * (p / (0.4 * rho) + 0.5 * pow(u, 2));;
  f[2] = u * (E + p);
}

FluxVar FluxVar::operator-(const FluxVar &rhs){
  return FluxVar(this->f[0] - rhs.f[0], this->f[1] - rhs.f[1], this->f[2] - rhs.f[2]);
}
