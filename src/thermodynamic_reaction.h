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

#ifndef BARITEREACTION_H
#define BARITEREACTION_H

#include "chemistry.h"
#include "uclib.h"
#include "math.h"
#include <cstdlib>
#include "equilibrium_formulation.h"
#include "activity_model.h"
#include "simple_reaction.h"
#include "simple_nucleation.h"

class ThermodynamicReaction
{
public:
    SimpleReaction &reaction;
    EquilibriumFormulation &equilibriumFormulation;
    ActivityModel &activityModel;
    SimpleNucleation &nucleationModel;

    // numerical
    const Real SMALL = 1e-16;

    // Constructor
    ThermodynamicReaction(SimpleReaction &reaction,
                          EquilibriumFormulation &equilibriumFormulation,
                          ActivityModel &activityModel,
                          SimpleNucleation &nucleationModel)
        : reaction(reaction),
          equilibriumFormulation(equilibriumFormulation),
          activityModel(activityModel),
          nucleationModel(nucleationModel)
    {
    }

    // Mean Molarity
    void MeanMolality(Real *meanMolality, int size, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            meanMolality[i] = MeanMolality(yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Equilibrium Function
    void Equilibrium(Real *equilibrium, int size, Real *Temperature)
    {
        for (int i = 0; i < size; i++)
        {
            equilibrium[i] = equilibriumFormulation.Equilibrium(Temperature[i]);
        }
    }

    // Ionic Strength
    void IonicStrength(Real *ionicStrength, int size, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            ionicStrength[i] = IonicStrength(yEtc1[i], yEtc2[i]);
        }
    }

    // Pitzers activitity coeff
    void ActivityCoefficient(Real *gamma, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            gamma[i] = activityModel.ActivityCoefficient(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Reaction rate of wall deposition in molality
    void WallReactionRateMolality(Real *wallReactionRate, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, Real *SaturationIndex, Real *WallDistance)
    {
        for (int i = 0; i < size; i++)
        {
            wallReactionRate[i] = nucleationModel.WallReactionRateMolality(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i], SaturationIndex[i], WallDistance[i]);
        }
    }

    // Reaction rate of wall deposition in mole fraction
    void WallReactionRateMoleFraction(Real *wallReactionRate, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, Real *SaturationIndex, Real *WallDistance)
    {
        Real mR;
        for (int i = 0; i < size; i++)
        {
            mR = nucleationModel.WallReactionRateMolality(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i], SaturationIndex[i], WallDistance[i]);
            wallReactionRate[i] = mR / TotalMolality(yEtc1[i], yEtc2[i]);
        }
    }

    // Nucleation Rates
    void NucleationRate(Real *nucleationRate, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, ThermodynamicReaction *obj, const Real (ThermodynamicReaction::*SaturationIndexFunction)(Real, Real, Real, Real, Real))
    {
        Real S;
        for (int i = 0; i < size; i++)
        {
            S = (obj->*SaturationIndexFunction)(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
            nucleationRate[i] = nucleationModel.NucleationRate(pow(S, 10), Temperature[i]);
        }
    }

    // Saturation Index Function
    void SaturationIndex(Real *saturationIndex, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            saturationIndex[i] = SaturationIndex(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Pitzer Wall Component A concentration
    void WallConcentrationA(Real *WallConcentration, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            WallConcentration[i] = WallConcentrationA(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Pitzer Wall Component B concentration
    void WallConcentrationB(Real *WallConcentration, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            WallConcentration[i] = WallConcentrationB(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Thermodynamic Equilibrium.
    void ActivityCorrectedEquilibrium(Real *correctedEquilibrium, int size, Real *equilibrium, Real *activity)
    {
        for (int i = 0; i < size; i++)
        {
            correctedEquilibrium[i] = pow(equilibrium[i] * activity[i], 1.0 / reaction.nu());
        }
    }

    // Returns Mean molality
    Real MeanMolality(Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        return MeanMolality(yA * mTot, yB * mTot);
    }

    Real TotalMolality(Real yEtc1, Real yEtc2)
    {
        return reaction.TotalMolality(yEtc1, yEtc2);
    }
    // Mean Concentration
    Real MeanMolality(Real mA, Real mB)
    {
        return reaction.MeanMolality(mA, mB);
    }

    // Iterative saturation to avoid overreaction - not in use
    Real SaturationIndex(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return log10(SaturationRate(T, yA, yB, yEtc1, yEtc2));
    }

    // Wall saturation model
    Real WallConcentrationA(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        double SR = SaturationRate(T, yA, yB, yEtc1, yEtc2);
        double mA = TotalMolality(yEtc1, yEtc2) * yA;
        return fmin(mA, 1 / SR * mA);
    }

    // Wall saturation model
    Real WallConcentrationB(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        double SR = SaturationRate(T, yA, yB, yEtc1, yEtc2);
        double mB = TotalMolality(yEtc1, yEtc2) * yB;
        return fmin(mB, 1 / SR * mB);
    }

    // Compute Ionic Strength
    Real IonicStrength(Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        const Real m[2] = {yEtc1 * mTot, yEtc2 * mTot};
        const Real Z[2] = {1, 2};
        return ChemistryFunctions::IonicStrength(m, Z, 2);
    }

    // Solve for saturation correcting for change in activity on reaction
    Real SaturationRate(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        Real mA = yA * mTot;
        Real mB = yB * mTot;

        if (fmin(mA, mB) < SMALL)
            return SMALL;

        const Real equilibrium = equilibriumFormulation.Equilibrium(T);
        const Real meanMolality = MeanMolality(mA, mB);
        const Real I = IonicStrength(yEtc1, yEtc2);
        const Real gamma = activityModel.ActivityCoefficient(T, yA, yB, yEtc1, yEtc2);

        return pow(meanMolality * gamma, reaction.nu()) / equilibrium;
    }

private:
};

#endif // BARITEREACTION_H
