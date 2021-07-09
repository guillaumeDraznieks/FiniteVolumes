#include "PrimitiveVar.h"
#include "FluxVar.h"
#include "appsolver.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include <cmath>

inline double kinetic_energy(double u)
{
    return 0.5 * pow(u, 2);
}

double AppSolver::solve(const PrimitiveVar &Wl, const PrimitiveVar &Wr, PrimitiveVar &Wsl, PrimitiveVar &Wsr, double zl, double zr){
	double xl , xr, alphal, alphar, prod, ps, us, taul,taur, el, er;
	xl = Wl.p/zl;
	xr = Wr.p/zr;
	alphal = zl/(zl+zr);
	alphar = 1-alphal;
	prod = zr*zl/(zr+zl);
	ps = prod*(xr+xl+Wl.u-Wr.u);
	us = alphal*Wl.u + alphar*Wr.u+ alphal*xl-alphar*xr;
	taul = Wl.tau+(us-Wl.u)/zl;
	taur = Wr.tau + (Wr.u-us)/zr;
	el = Wl.e + kinetic_energy(Wl.u)-kinetic_energy(us)+(Wl.p*Wl.u-ps*us)/zl; // Energie interne, pas totale!!
	er = Wr.e + kinetic_energy(Wr.u)-kinetic_energy(us)+(ps*us-Wl.p*Wl.u)/zl; // Energie interne, pas totale!!
	Wsl.lag_set(taul, us,el*PrimitiveVar::G8/taul);
	Wsr.lag_set(taur, us,er*PrimitiveVar::G8/taur);
	return ps;

}
