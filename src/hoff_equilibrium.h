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

#ifndef HOFF_EQUILIBRIUM_H
#define HOFF_EQUILIBRIUM_H

#include "chemistry.h"
#include "uclib.h"
#include "math.h"
#include <cstdlib>
#include "equilibrium_formulation.h"

class HoffEquilibrium : public EquilibriumFormulation
{
public:
    const Real log_k; // FROM PREEQC
    const Real delta_h;
    const Real T0;

    HoffEquilibrium(Real log_k, Real delta_h, Real T0) : log_k(log_k), delta_h(delta_h), T0(T0)
    {
        // const Real log_k = -9.87; // FROM PREEQC
        // const Real delta_h = 6.35 * 4186.80;
        // const Real T0 = ChemistryFunctions::T0()
    }

    // Equilibrium concentration at T
    Real Equilibrium(Real T)
    {
        return pow(10.0, log_k + delta_h / ChemistryFunctions::R() * (1 / T0 - 1 / T));
    }
};

#endif // HOFF_EQUILIBRIUM_H
