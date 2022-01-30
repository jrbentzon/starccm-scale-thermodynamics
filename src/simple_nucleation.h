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

#ifndef SIMPLE_NUCLEATION_H
#define SIMPLE_NUCLEATION_H

#include "chemistry.h"
#include "uclib.h"
#include "math.h"
#include <cstdlib>
#include "simple_reaction.h"

class SimpleNucleation
{
public:
    SimpleReaction *reaction;

    // Nucleation
    const Real MolarMass;
    const Real sigma;
    const Real theta;
    const Real rho;
    const Real A_n;

    // numerical
    const Real SMALL = 1e-16;
    const Real pi = M_PI;
    const Real DiffusionLayerThickness = 1e-3; // Random number - 1mm need to be linked to diffusion rate
    const Real TransitionLayerThickness = 2e-4;
    const Real DrivingMolalitySmoothingGap = 1e-8; // molality used for derivatives - not a physical property

    SimpleNucleation(
        SimpleReaction *reaction,
        Real MolarMass,
        Real sigma,
        Real theta,
        Real rho,
        Real A_n)
        : MolarMass(MolarMass),
          sigma(sigma),
          theta(theta),
          rho(rho),
          A_n(A_n)
    {
        this->reaction = reaction;
    }

    // Contact Angle Function
    const Real f()
    {
        return (2.0 - 3.0 * cos(theta) + pow(cos(theta), 3.0)) / 8.0;
    }

    // Nucleation rate
    const Real NucleationRate(Real S, Real T)
    {
        if (S <= 1.001)
        {
            return 0;
        }

        const Real Vm = MolarMass / rho;
        const Real ExponentNumerator = 16.0 * ChemistryFunctions::N_A() * pi * pow(Vm, 2.0) * pow(sigma, 3.0);
        const Real ExponentDenominator = 3.0 * pow(ChemistryFunctions::R() * T, 3.0) * pow(log(S), 2.0);

        return A_n * exp(-ExponentNumerator / ExponentDenominator);
    }

    // Nucleation Critical Radius
    const Real CriticalRadius(Real S, Real T)
    {
        if (S <= 1.001)
        {
            return 0;
        }

        const Real ExponentNumerator = 2.0 * MolarMass * sigma;
        const Real ExponentDenominator = rho * ChemistryFunctions::R() * T * log(S);

        return exp(-ExponentNumerator / ExponentDenominator);
    }

    // Produced molarity temporal derivative
    const Real BulkReactionRate(Real S, Real T, Real WallDistance) // mass per time
    {
        Real dNdt = NucleationRate(S, T);
        Real dr = CriticalRadius(S, T);
        Real dV = 4.0 / 3.0 * pi * pow(dr, 3);

        return (1.0 - wallReactionRatio(WallDistance)) * dV * rho * dNdt;
    }

    // How much of this cells content can react on the wall \in [0,1]
    const Real wallReactionRatio(Real WallDistance)
    {
        if (WallDistance < DiffusionLayerThickness - 0.5 * TransitionLayerThickness)
        {
            return 1;
        }
        else if (WallDistance < DiffusionLayerThickness + 0.5 * TransitionLayerThickness)
        {
            return 1 - (WallDistance - (DiffusionLayerThickness - 0.5 * TransitionLayerThickness)) / TransitionLayerThickness;
        }
        else
        {
            return 0;
        }
    }

    // Per cell functions:
    // Molality Based Wall Reaction Rate
    const Real WallReactionRateMolality(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2, Real SaturationIndex, Real WallDistance)
    {
        const Real mTot = reaction->TotalMolality(yEtc1, yEtc2);
        const Real mA = yA * mTot;
        const Real mB = yB * mTot;
        const Real meanMolality = reaction->MeanMolality(mA, mB);
        const Real SR = pow(SaturationIndex, 10);
        if (SR - 1.0 < SMALL)
            return 0;
        else
            return WallReactionRateMolality(SR, mA, mB, meanMolality, WallDistance);
    }

    // Molality Based Wall Reaction Rate
    const Real WallReactionRateMolality(Real SaturationRate, Real mA, Real mB, Real meanMolality, Real WallDistance)
    {
        const Real kWall = wallReactionRatio(WallDistance); // int_0,yScale y*S dy
        return kWall * (SaturationRate - 1) / SaturationRate * DrivingMolality(mA, mB, meanMolality);
    }

    // Driving molality - index
    const int DrivingMolalityIndex(Real mA, Real mB, Real mMean)
    {
        if (mA < mB && mA < mMean)
            return 0;
        else if (mB < mMean)
            return 1;
        else
            return 2;
    }

    // Driving molality in molality
    const Real DrivingMolality(Real mA, Real mB, Real mMean)
    {
        const int index = DrivingMolalityIndex(mA, mB, mMean);
        if (index == 0)
            return mA;
        else if (index == 1)
            return mB;
        else
            return mMean;
    }

    // Excess of a given specie in reaction
    const Real ExcessReactionMolality(Real m_i, Real drivingMolality)
    {
        return m_i - drivingMolality;
    }

    // Rate of how driving the molality 0 if not dribing, 1 if driving and smooth in gap
    const Real DrivingMolalityRate(Real excessMolality)
    {
        return smoothstep(DrivingMolalitySmoothingGap, 0, excessMolality);
    }

    const Real DrivingMolalityRate(Real m_i, Real drivingMolality)
    {
        return DrivingMolalityRate(ExcessReactionMolality(m_i, drivingMolality));
    }

    // smoothstep function 0 at edge 0, 1 at edge 1 - polynomial in between - smooth to 1st order
    Real smoothstep(Real edge0, Real edge1, Real x)
    {
        // Scale, bias and saturate x to 0..1 range
        x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
        // Evaluate polynomial
        return x * x * (3 - 2 * x);
    }

    Real clamp(Real x, Real lowerlimit, Real upperlimit)
    {
        if (x < lowerlimit)
            x = lowerlimit;
        if (x > upperlimit)
            x = upperlimit;
        return x;
    }

    // Derivatives:

    // Wall Reaction to Saturation
    const Real dRw_dSR(Real S, Real drivingMolality)
    {
        // d(S-1)/dS + d(S^-1)/dS
        return (1.0 - (1.0 - S) * pow(S, -2)) * drivingMolality;
    }

    // Saturation to molality
    const Real dSR_dmMean(Real gamma, Real K)
    {
        return gamma / K;
    }

    const Real dmMean_dmi(Real mMean, Real m_i, Real nu_i)
    {
        if (m_i < SMALL)
            return 0;

        return mMean * nu_i / (reaction->nu() * m_i);
    }

    // SaturationRate to activity
    const Real dSR_dgamma(Real mMean, Real K)
    {
        // TODO: find derivative
        return mMean / K;
    }

    //ionic strength to molality
    const Real dI_dmi(Real Z_i)
    {
        return 0.5 * Z_i * Z_i;
    }

    // Molality to mole fraction
    const Real dm_dyi(Real y_i, Real mTot)
    {
        if (y_i < SMALL)
            return 0;

        const Real m_i = y_i * mTot;
        return pow(y_i - 1, -2) * (mTot - m_i);
    }

    // Chain rule:
    // Ionic strength to mole fraciton
    const Real dI_dyi(Real y_i, Real Z_i, Real mTot)
    {
        return dI_dmi(Z_i) * dm_dyi(y_i, mTot);
    }
};

#endif // SIMPLE_NUCLEATION_H
