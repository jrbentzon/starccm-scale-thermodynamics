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
#include "thermodynamic_reaction.h"
#include "simple_reaction.h"
#include "activity_model.h"
#include "pitzer_activity_model.h"
#include "simple_nucleation.h"

using namespace std;

const ThermodynamicReaction getReaction()
{
    // Reaction Model
    const Real nu_A = 1;
    const Real nu_B = 1;
    const Real nu_P = -1;
    const Real Z_A = 2;
    const Real Z_B = -2;

    const SimpleReaction reactionModel = SimpleReaction(nu_A, nu_B, nu_P, Z_A, Z_B);

    // Equilibrium Model
    const Real log_k = -9.87;
    const Real delta_h = 6.35 * 4186.80;
    const Real T0 = ChemistryFunctions::T0();
    const EquilibriumFormulation equilibriumModel = HoffEquilibrium(log_k, delta_h, T0);

    // Acticity Coefficient Model
    const Real beta_0 = 0;
    const Real beta_1 = 0;
    const Real beta_2 = 0;
    const Real C_phi = 0;
    const ActivityModel activityModel = PitzerActivityModel(reactionModel, beta_0, beta_1, beta_2, C_phi);

    // Nucleation Model
    const Real MolarMass = 0.1000894;
    const Real sigma = 40e-3;
    const Real theta = 10.0 / 180.0 * M_PI;
    const Real rho = 2710;
    const Real A_n = 1;
    SimpleNucleation nucleationModel = SimpleNucleation(reactionModel, MolarMass, sigma, theta, rho, A_n);

    return ThermodynamicReaction(reactionModel, equilibriumModel, activityModel, nucleationModel);
}

void EquilibriumConstant(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    ThermodynamicReaction React = getReaction();
    React.Equilibrium(result, size, Temperature);
}

void PitzerActivity(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    ThermodynamicReaction React = getReaction();
    React.ActivityCoefficient(result, size, Temperature, yA, yB, yEtc_1, yEtc_2);
}
void PitzerSaturationIndex(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    ThermodynamicReaction React = getReaction();
    React.SaturationIndex(result, size, Temperature, yA, yB, yEtc_1, yEtc_2);
}

void IonicStrength(Real *result, int size, Real *yEtc_1, Real *yEtc_2)
{
    ThermodynamicReaction React = getReaction();
    React.IonicStrength(result, size, yEtc_1, yEtc_2);
}

void MeanMolality(Real *result, int size, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    ThermodynamicReaction React = getReaction();
    React.MeanMolality(result, size, yA, yB, yEtc_1, yEtc_2);
}

void uclib()
{
    ucfunc((void *)EquilibriumConstant, "ScalarFieldFunction", "Equilibrium Constant");
    ucarg((void *)EquilibriumConstant, "Cell", "Temperature", sizeof(Real));

    ucfunc((void *)IonicStrength, "ScalarFieldFunction", "Ionic Strength");
    ucarg((void *)IonicStrength, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)IonicStrength, "Cell", "$yEtc_2-", sizeof(Real));

    ucfunc((void *)MeanMolality, "ScalarFieldFunction", "Mean Molality");
    ucarg((void *)MeanMolality, "Cell", "$yBa_2+", sizeof(double));
    ucarg((void *)MeanMolality, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)MeanMolality, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)MeanMolality, "Cell", "$yEtc_2-", sizeof(Real));

    ucfunc((void *)PitzerActivity, "ScalarFieldFunction", "Pitzer Activity Coefficient");
    ucarg((void *)PitzerActivity, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)PitzerActivity, "Cell", "$yBa_2+", sizeof(double));
    ucarg((void *)PitzerActivity, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)PitzerActivity, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)PitzerActivity, "Cell", "$yEtc_2-", sizeof(Real));

    ucfunc((void *)PitzerSaturationIndex, "ScalarFieldFunction", "Pitzer Saturation Index");
    ucarg((void *)PitzerSaturationIndex, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)PitzerSaturationIndex, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)PitzerSaturationIndex, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)PitzerSaturationIndex, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)PitzerSaturationIndex, "Cell", "$yEtc_2-", sizeof(Real));
}
