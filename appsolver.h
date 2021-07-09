#ifndef APPSOLVER_H
#define APPSOLVER_H

#include <vector>

class AppSolver
{
public:

    static double solve(const PrimitiveVar &Wl, const PrimitiveVar &Wr, PrimitiveVar &Wsl, PrimitiveVar &Wsr, double zl, double zr);


};

#endif // APPSOLVER_H