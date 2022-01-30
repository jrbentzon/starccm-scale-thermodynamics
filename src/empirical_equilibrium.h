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

#ifndef EMPIRICAL_EQUILIBRIUM_H
#define EMPIRICAL_EQUILIBRIUM_H

#include "chemistry.h"
#include "uclib.h"
#include "math.h"
#include <cstdlib>
#include "equilibrium_formulation.h"

class EmpiricalEquilibrium : public EquilibriumFormulation
{
public:
    const Real analytical_expression[4];

    EmpiricalEquilibrium(Real A, Real B, Real C, Real D) : analytical_expression({A, B, C, D})
    {
    }

    // Equilibrium concentration at T
    Real Equilibrium(Real T)
    {
        return pow(10.0, analytical_expression[0] + analytical_expression[1] * T + analytical_expression[2] / T + analytical_expression[3] * log10(T));
    }
};

#endif // EMPIRICAL_EQUILIBRIUM_H
