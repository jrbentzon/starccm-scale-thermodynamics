/*
MIT License

Copyright (c) 2022 Jakob Roar Bentzon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef CHEMFUNC_H
#define CHEMFUNC_H

namespace ChemistryFunctions
{

static Real k_b()
{
    return 1.38064852e-23;
}

static Real N_A()
{
    return 6.022140857e23;
}

static Real MolarMassOfWater()
{
    return 18.01528e-3;
}

static Real R()
{
    return N_A() * k_b(); // should be precompiled
}

static Real T0()
{
    return 273.15 + 25;
}

static Real electronicCharge()
{
    return 1.60206e-19;
}

static Real permittivityVacuum()
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
