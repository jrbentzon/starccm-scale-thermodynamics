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

#ifndef SIMPLE_REACTION_H
#define SIMPLE_REACTION_H

#include "chemistry.h"
#include "uclib.h"
#include "math.h"
#include <cstdlib>

class SimpleReaction
{
public:
    Real nu_A;
    Real nu_B;
    Real nu_P;
    Real Z_A;
    Real Z_B;

    Real SMALL = 1e-16;

    SimpleReaction(Real nu_A, Real nu_B, Real nu_P, Real Z_A, Real Z_B)
        : nu_A(nu_A),
          nu_B(nu_B),
          nu_P(nu_P),
          Z_A(Z_A),
          Z_B(Z_B)
    {
    }

    Real nu()
    {
        return nu_A + nu_B;
    }

    // Total Concentration
    Real TotalMolality(Real yEtc1, Real yEtc2)
    {
        return pow(ChemistryFunctions::MolarMassOfWater() / (1 - (yEtc1 + yEtc2)) + SMALL, -1);
    }

    // Returns Mean molality
    Real MeanMolality(Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        return MeanMolality(yA * mTot, yB * mTot);
    }

    // Mean Concentration
    Real MeanMolality(Real mA, Real mB)
    {
        return pow(pow(fmax(mA, SMALL), nu_A) * pow(fmax(mB, SMALL), nu_B) + SMALL, 1.0 / nu());
    }
};

#endif // SIMPLE_REACTION_H
