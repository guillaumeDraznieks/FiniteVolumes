#ifndef FLUXVAR_H
#define FLUXVAR_H


class FluxVar
{
  public:
    double f[3];

    FluxVar(double a, double b, double c);
    FluxVar();

    void set (double density, double velocity, double pressure);
    void lag_set (double u, double p);

    FluxVar operator-(const FluxVar &rhs);
    //~FluxVar();


};

#endif //FLUXVAR_H