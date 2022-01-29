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

#include "barite.H"

using namespace std;

const BariteReaction getBaSO4Reaction()
{
    const Real nu_A = 1;
    const Real nu_B = 1;
    const Real nu_P = -1;
    const Real Z_A = 2;
    const Real Z_B = -2;
    return BariteReaction(nu_A, nu_B, nu_P, Z_A, Z_B);
}

void DebyeHuckelActivity(Real *result, int size, Real *Temperature, Real *mA, Real *mB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.DebyeHuckelActivity(result, size, Temperature, yEtc_1, yEtc_2);
}

void DebyeHuckelSaturationIndex(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.DebyeHuckelSaturationIndex(result, size, Temperature, yA, yB, yEtc_1, yEtc_2);
}

void EquilibriumConstant(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.Equilibrium(result, size, Temperature);
}

void PitzerActivity(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.pitzerActivityCoefficient(result, size, Temperature, yA, yB, yEtc_1, yEtc_2);
}
void PitzerSaturationIndex(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.PitzerSaturationIndex(result, size, Temperature, yA, yB, yEtc_1, yEtc_2);
}

void IonicStrength(Real *result, int size, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.IonicStrength(result, size, yEtc_1, yEtc_2);
}

void MeanMolality(Real *result, int size, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.MeanMolality(result, size, yA, yB, yEtc_1, yEtc_2);
}

void ExtendedDebyeHuckelNucleationRate(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    const Real (NucleationRate::*SatFunctionPointer)(Real, Real, Real, Real, Real) = &BariteReaction::ExtendedDebyeHuckelSaturationIndex;
    React.NucleationRate(result, size, Temperature, yA, yB, yEtc_1, yEtc_2, &React, SatFunctionPointer);
}

void ExtendedDebyeHuckelWallReactionRateMolality(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2, Real *SaturationIndex, Real *WallDistance)
{
    BariteReaction React = getBaSO4Reaction();
    React.WallReactionRateMolality(result, size, Temperature, yA, yB, yEtc_1, yEtc_2, SaturationIndex, WallDistance);
}

void ExtendedDebyeHuckelWallReactionRateMoleFraction(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2, Real *SaturationIndex, Real *WallDistance)
{
    BariteReaction React = getBaSO4Reaction();
    React.WallReactionRateMoleFraction(result, size, Temperature, yA, yB, yEtc_1, yEtc_2, SaturationIndex, WallDistance);
}

void ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2, Real *SaturationIndex, Real *WallDistance)
{
    BariteReaction React = getBaSO4Reaction();
    React.DebyeHuckelWallReactionDerivativeA(result, size, Temperature, yA, yB, yEtc_1, yEtc_2, SaturationIndex, WallDistance);
}

void ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2, Real *SaturationIndex, Real *WallDistance)
{
    BariteReaction React = getBaSO4Reaction();
    React.DebyeHuckelWallReactionDerivativeB(result, size, Temperature, yA, yB, yEtc_1, yEtc_2, SaturationIndex, WallDistance);
}

void WallConcentrationA(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.PitzerWallConcentrationA(result, size, Temperature, yA, yB, yEtc_1, yEtc_2);
}

void WallConcentrationB(Real *result, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc_1, Real *yEtc_2)
{
    BariteReaction React = getBaSO4Reaction();
    React.PitzerWallConcentrationB(result, size, Temperature, yA, yB, yEtc_1, yEtc_2);
}

void uclib()
{

    ucfunc((void *)EquilibriumConstant, "ScalarFieldFunction", "Equilibrium Constant");
    ucarg((void *)EquilibriumConstant, "Cell", "Temperature", sizeof(Real));

    ucfunc((void *)DebyeHuckelActivity, "ScalarFieldFunction", "Debye Huckel Activity Coefficient");
    ucarg((void *)DebyeHuckelActivity, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)DebyeHuckelActivity, "Cell", "$yBa_2+", sizeof(double));
    ucarg((void *)DebyeHuckelActivity, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)DebyeHuckelActivity, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)DebyeHuckelActivity, "Cell", "$yEtc_2-", sizeof(Real));

    ucfunc((void *)DebyeHuckelSaturationIndex, "ScalarFieldFunction", "Debye Huckel Saturation Index");
    ucarg((void *)DebyeHuckelSaturationIndex, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)DebyeHuckelSaturationIndex, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)DebyeHuckelSaturationIndex, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)DebyeHuckelSaturationIndex, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)DebyeHuckelSaturationIndex, "Cell", "$yEtc_2-", sizeof(Real));

    ucfunc((void *)ExtendedDebyeHuckelWallReactionRateMolality, "ScalarFieldFunction", "Debye Huckel Wall Deposition (Molality)");
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMolality, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMolality, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMolality, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMolality, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMolality, "Cell", "$yEtc_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMolality, "Cell", "$UserDebyeHuckelSaturationIndex", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMolality, "Cell", "$WallDistance", sizeof(Real));

    ucfunc((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "ScalarFieldFunction", "Debye Huckel Wall Deposition (Mole Fraction)");
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "Cell", "$yEtc_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "Cell", "$UserDebyeHuckelSaturationIndex", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateMoleFraction, "Cell", "$WallDistance", sizeof(Real));

    ucfunc((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "ScalarFieldFunction", "Debye Huckel Wall Deposition Derivative Ba (Mole Fraction)");
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "Cell", "$yEtc_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "Cell", "$UserDebyeHuckelSaturationIndex", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeAMoleFraction, "Cell", "$WallDistance", sizeof(Real));

    ucfunc((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "ScalarFieldFunction", "Debye Huckel Wall Deposition Derivative SO4 (Mole Fraction)");
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "Cell", "$yEtc_2-", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "Cell", "$UserDebyeHuckelSaturationIndex", sizeof(Real));
    ucarg((void *)ExtendedDebyeHuckelWallReactionRateDerivativeBMoleFraction, "Cell", "$WallDistance", sizeof(Real));

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

    ucfunc((void *)WallConcentrationA, "ScalarFieldFunction", "Pitzer Diff Model Wall mBa_2+");
    ucarg((void *)WallConcentrationA, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)WallConcentrationA, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)WallConcentrationA, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)WallConcentrationA, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)WallConcentrationA, "Cell", "$yEtc_2-", sizeof(Real));

    ucfunc((void *)WallConcentrationB, "ScalarFieldFunction", "Pitzer Diff Model Wall mSO4_2-");
    ucarg((void *)WallConcentrationB, "Cell", "Temperature", sizeof(Real));
    ucarg((void *)WallConcentrationB, "Cell", "$yBa_2+", sizeof(Real));
    ucarg((void *)WallConcentrationB, "Cell", "$ySO4_2-", sizeof(Real));
    ucarg((void *)WallConcentrationB, "Cell", "$yEtc_1-", sizeof(Real));
    ucarg((void *)WallConcentrationB, "Cell", "$yEtc_2-", sizeof(Real));
}
