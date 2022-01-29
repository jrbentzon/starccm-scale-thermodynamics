#ifndef CHEMFUNC_H
#define CHEMFUNC_H

// #include "uclib.h"
// #include "math.h"
// #include "stdlib.h"

namespace ChemistryFunctions
{

static const Real k_b()
{

    return 1.38064852e-23;
}

static const Real N_A()
{
    return 6.022140857e23;
}

static const Real MolarMassOfWater()
{
    return 18.01528e-3;
}

static const Real R()
{
    return N_A() * k_b(); // should be precompiled
}

static const Real T0()
{
    return 273.15 + 25;
}

static const Real electronicCharge()
{
    return 1.60206e-19;
}

static const Real permittivityVacuum()
{
    return 8.8542e-12;
}

static Real permittivityWater()
{
    //Make function of T
    return 78.4;
}

static Real densityWater()
{
    //Make function of T
    return 997;
}

static Real IonicStrength(const Real *m, const Real *Z, const int noReactants)
{
    Real I = 0;
    for (int j = 0; j < noReactants; j++)
    {
        I = I + m[j] * pow(Z[j], 2); //0.5 should be taken outside loop for better performance
    }

    I = 0.5 * I;
    // AVOID zero
    if (I < 1e-18)
    {
        I = 1e-18;
    }
    return I;
}

static Real unitMolar()
{
    return 1;
}

}; // namespace ChemistryFunctions

#endif // CHEMFUNC_H
