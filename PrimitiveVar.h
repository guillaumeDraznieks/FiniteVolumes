#ifndef PRIMITIVEVAR_H
#define PRIMITIVEVAR_H

#include <vector>

class PrimitiveVar
{
  public:

    static double TOL,G0,G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12;

    double tau, rho, u, p;
    double a, e, a_lag;

    PrimitiveVar(double density, double velocity, double pressure);
    PrimitiveVar();
    PrimitiveVar(const PrimitiveVar &W);

    void set(double density, double velocity, double pressure);
    void lag_set(double tau, double velocity, double pressure);
    void set(const PrimitiveVar &W);

    static double E(double density, double velocity, double pressure);
    static double E(const PrimitiveVar &W);
    static double kinetic_energy(double u);
    static double sound_speed(double p, double rho);
    static double lag_sound_speed(double p, double tau);
    static double internal_energy(double p, double rho);
};

#endif //PRIMITIVEVAR_H