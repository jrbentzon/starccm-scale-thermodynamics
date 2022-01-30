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

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include "uclib.h"
#include "equilibrium_formulation.h"
#include "hoff_equilibrium.h"
#include "chemistry.h"
#include "simple_reaction.h"
#include "activity_model.h"
#include "pitzer_activity_model.h"

using namespace std;

void EquilibriumConstant(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    // Equilibrium Model
    const Real log_k = -9.87;
    const Real delta_h = 6.35 * 4186.80;
    const Real T0 = ChemistryFunctions::T0();
    HoffEquilibrium equilibriumModel = HoffEquilibrium(log_k, delta_h, T0);

    for (int i = 0; i < size; i++)
    {
        result[i] = equilibriumModel.Equilibrium(Temperature[i]);
    }
}

void PitzerActivity(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    // Acticity Coefficient Model
    static Real nu_A = 1;
    static Real nu_B = 1;
    static Real nu_P = -1;
    static Real Z_A = 2;
    static Real Z_B = -2;
    static Real beta_0 = 0;
    static Real beta_1 = 0;
    static Real beta_2 = 0;
    static Real C_phi = 0;
    PitzerActivityModel activityModel = PitzerActivityModel(nu_A, nu_B, Z_A, Z_B, beta_0, beta_1, beta_2, C_phi);
    for (int i = 0; i < size; i++)
    {
        result[i] = activityModel.ActivityCoefficient(Temperature[i], yA[i], yB[i], yEtc_1[i], yEtc_2[i]);
    }
}

void uclib()
{
    ucfunc((void *)EquilibriumConstant, "ScalarFieldFunction", "Equilibrium Constant");
    ucarg((void *)EquilibriumConstant, "Cell", "Temperature", sizeof(Real));

    ucfunc((void *)PitzerActivity, "ScalarFieldFunction", "Pitzer Activity Coefficient");
    ucarg((void *)PitzerActivity, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)PitzerActivity, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)PitzerActivity, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)PitzerActivity, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)PitzerActivity, "Cell", "$yEtc_2-", sizeof(Real));
}
