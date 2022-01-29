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

#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include "functions.H"
#include "uclib.h"
#include "math.h"
#include <cstdlib>

class BariteReaction
{
public:
    const Real sigma;
    const Real theta;
    const Real rho;
    const Real An;

    const Real nu_A;
    const Real nu_B;
    const Real nu_P;
    const Real nu;

    const Real Z_A, Z_B;

    const Real log_k = -9.87; // FROM PREEQC
    const Real delta_h = 6.35 * 4186.80;
    const Real analytical_expression[4] = {-282.43, -8.972e-2, 5822, 113.08};

    const Real MolarMass;

    const Real beta_0;
    const Real beta_1;
    const Real beta_2;
    const Real C_Phi;

    const Real T0;
    const Real R;
    const Real SMALL;
    const Real pi;
    const Real rho_w;
    const Real e;
    const Real eps_0;
    const Real eps_r;
    const Real N_A;
    const Real k_b;

    const Real DiffusionLayerThickness = 1e-3; // Random number - 1mm need to be linked to diffusion rate
    const Real TransitionLayerThickness = 2e-4;
    const Real DrivingMolalitySmoothingGap = 1e-8; // molality used for derivatives - not a physical property

    // Constructor
    BariteReaction(Real nu_A, Real nu_B, Real nu_P, Real Z_A, Real Z_B) : SMALL(1e-16),
                                                                    R(ChemistryFunctions::R()),
                                                                    T0(ChemistryFunctions::T0()),
                                                                    pi(3.14159265358979323846),
                                                                    N_A(ChemistryFunctions::N_A()),
                                                                    e(ChemistryFunctions::electronicCharge()),
                                                                    rho_w(ChemistryFunctions::densityWater()),
                                                                    eps_0(ChemistryFunctions::permittivityVacuum()),
                                                                    eps_r(ChemistryFunctions::permittivityWater()),
                                                                    k_b(ChemistryFunctions::k_b()),
                                                                    nu_A(nu_A),
                                                                    nu_B(nu_B),
                                                                    nu_P(nu_P),
                                                                    nu(nu_A + nu_B),
                                                                    Z_A(Z_A),
                                                                    Z_B(Z_B),
                                                                    beta_0(0),
                                                                    beta_1(0),
                                                                    beta_2(0),
                                                                    C_Phi(0),
                                                                    MolarMass(0.1000894),
                                                                    sigma(40e-3),
                                                                    theta(10.0 / 180.0 * pi),
                                                                    rho(2710),
                                                                    An(1)
    {
    }

    // Mean Molarity
    const void MeanMolality(Real *meanMolality, int size, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            meanMolality[i] = MeanMolality(yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Contact Angle Function
    const Real f()
    {
        return (2.0 - 3.0 * cos(theta) + pow(cos(theta), 3.0)) / 8.0;
    }

    // Equilibrium Function
    const void Equilibrium(Real *equilibrium, int size, Real *Temperature)
    {
        for (int i = 0; i < size; i++)
        {
            equilibrium[i] = Equilibrium(Temperature[i]);
        }
    }

    // Ionic Strength
    const void IonicStrength(Real *ionicStrength, int size, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            ionicStrength[i] = IonicStrength(yEtc1[i], yEtc2[i]);
        }
    }

    // Debye Huckel Derivative of Wall Reaction to mA
    const void DebyeHuckelWallReactionDerivativeA(Real *wallReactionRateDerivativeA, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, Real *SaturationIndex, Real *WallDistance)
    {
        for (int i = 0; i < size; i++)
        {
            wallReactionRateDerivativeA[i] = DebyeHuckelWallReactionDerivativeA(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i], WallDistance[i], SaturationIndex[i]);
        }
    }

    // Debye Huckel Derivative of Wall Reaction to mB
    const void DebyeHuckelWallReactionDerivativeB(Real *wallReactionRateDerivativeA, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, Real *SaturationIndex, Real *WallDistance)
    {
        for (int i = 0; i < size; i++)
        {
            wallReactionRateDerivativeA[i] = DebyeHuckelWallReactionDerivativeA(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i], WallDistance[i], SaturationIndex[i]);
        }
    }

    // Pitzers activitity coeff
    const void pitzerActivityCoefficient(Real *gamma, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            gamma[i] = pitzerActivityCoefficient(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Reaction rate of wall deposition in molality
    const void WallReactionRateMolality(Real *wallReactionRate, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, Real *SaturationIndex, Real *WallDistance)
    {
        for (int i = 0; i < size; i++)
        {
            wallReactionRate[i] = WallReactionRateMolality(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i], SaturationIndex[i], WallDistance[i]);
        }
    }

    // Reaction rate of wall deposition in mole fraction
    const void WallReactionRateMoleFraction(Real *wallReactionRate, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, Real *SaturationIndex, Real *WallDistance)
    {
        Real mR;
        for (int i = 0; i < size; i++)
        {
            mR = WallReactionRateMolality(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i], SaturationIndex[i], WallDistance[i]);
            wallReactionRate[i] = mR / TotalMolality(yEtc1[i], yEtc2[i]);
        }
    }

    // Nucleation Rates
    const void NucleationRate(Real *nucleationRate, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2, reaction *obj, const Real (BariteReaction::*SaturationIndexFunction)(Real, Real, Real, Real, Real))
    {
        Real S;
        for (int i = 0; i < size; i++)
        {
            S = (obj->*SaturationIndexFunction)(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
            nucleationRate[i] = NucleationRate(pow(S, 10), Temperature[i]);
        }
    }

    // Saturation Index Function
    const void PitzerSaturationIndex(Real *saturationIndex, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            saturationIndex[i] = IterativePitzerSaturationIndex(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

     // Pitzer Wall Component A concentration
    const void PitzerWallConcentrationA(Real *WallConcentrationA, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            WallConcentrationA[i] = PitzerWallConcentrationA(Temperature[i],yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

     // Pitzer Wall Component A concentration
    const void PitzerWallConcentrationB(Real *WallConcentrationA, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            WallConcentrationA[i] = PitzerWallConcentrationB(Temperature[i],yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Saturation Index Function
    const void DebyeHuckelSaturationIndex(Real *saturationIndex, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            saturationIndex[i] = DebyeHuckelSaturationIndex(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Saturation Index Function
    const void ExtendedDebyeHuckelSaturationIndex(Real *saturationIndex, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            saturationIndex[i] = ExtendedDebyeHuckelSaturationIndex(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Thermodynamic Equilibrium.
    const void ActivityCorrectedEquilibrium(Real *correctedEquilibrium, int size, Real *equilibrium, Real *activity)
    {

        for (int i = 0; i < size; i++)
        {
            correctedEquilibrium[i] = pow(equilibrium[i] * activity[i], 1.0 / nu);
        }
    }

    // Extended Debye Huckel Activity Coeff
    const void ExtendedDebyeHuckelActivity(Real *DebyeHuckelAct, int size, Real *Temperature, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            DebyeHuckelAct[i] = ExtendedDebyeHuckelActivity(Temperature[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Standard Debye Huckel Activity Coeff
    const void DebyeHuckelActivity(Real *DebyeHuckelAct, int size, Real *Temperature, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            DebyeHuckelAct[i] = DebyeHuckelActivity(Temperature[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Per cell functions:

    // Molality Based Wall Reaction Rate
    const Real WallReactionRateMolality(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2, Real SaturationIndex, Real WallDistance)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        const Real mA = yA * mTot;
        const Real mB = yB * mTot;
        const Real meanMolality = MeanMolality(mA, mB);
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

    // Derivative of wall reaction rate to reactant A
    const Real DebyeHuckelWallReactionDerivativeA(Real Temperature, Real y_A, Real y_B, Real yEtc1, Real yEtc2, Real WallDistance, Real SaturationIndex)
    {
        const Real kwr = wallReactionRatio(WallDistance);
        if (kwr > SMALL) // Save computation for non-wall cells
        {

            const Real mMean = MeanMolality(y_A, y_B, yEtc1, yEtc2);
            const Real mTot = TotalMolality(yEtc1, yEtc2);
            const Real m_A = mTot * y_A, m_B = mTot * y_B;
            const Real K = Equilibrium(Temperature);
            const Real gamma = ExtendedDebyeHuckelActivity(Temperature, yEtc1, yEtc2);
            const Real drivingMolality = DrivingMolality(m_A, m_B, mMean);
            const Real SR = pow(SaturationIndex, 10);

            return kwr * dRw_dyi(Temperature, y_A, m_A, nu_A, Z_A, mMean, SR, K, gamma, mTot, drivingMolality);
        }
        else
            return 0;
    }

    // Derivative of wall reaction rate to reactant B
    const Real DebyeHuckelWallReactionDerivativeB(Real Temperature, Real y_A, Real y_B, Real yEtc1, Real yEtc2, Real WallDistance, Real SaturationIndex)
    {
        const Real kwr = wallReactionRatio(WallDistance);
        if (kwr > SMALL || y_A < SMALL) // Save computation for non-wall cells
        {

            const Real mMean = MeanMolality(y_A, y_B, yEtc1, yEtc2);
            const Real mTot = TotalMolality(yEtc1, yEtc2);
            const Real m_A = mTot * y_A, m_B = mTot * y_B;
            const Real K = Equilibrium(Temperature);
            const Real gamma = ExtendedDebyeHuckelActivity(Temperature, yEtc1, yEtc2);
            const Real drivingMolality = DrivingMolality(m_A, m_B, mMean);
            const Real SR = pow(SaturationIndex, 10);

            return kwr * dRw_dyi(Temperature, y_B, m_B, nu_B, Z_B, mMean, SR, K, gamma, mTot, drivingMolality);
        }
        else
            return 0;
    }

    // Returns Mean molality
    const Real MeanMolality(Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        return MeanMolality(yA * mTot, yB * mTot);
    }

    const Real TotalMolality(Real yEtc1, Real yEtc2)
    {
        return pow(ChemistryFunctions::MolarMassOfWater() / (1 - (yEtc1 + yEtc2)) + SMALL,-1);
    }
    // Mean Concentration
    const Real MeanMolality(Real mA, Real mB)
    {
        return pow(pow(fmax(mA, SMALL), nu_A) * pow(fmax(mB, SMALL), nu_B) + SMALL, 1.0 / nu);
    }


 // Equilibrium concentration at T
    const Real Equilibrium(Real T)
    {
        return EmpiricalEquilibrium(T);
    }

    // Equilibrium concentration at T
    const Real HoffEquilibrium(Real T)
    {
        return pow(10.0, log_k + delta_h / R * (1 / T0 - 1 / T));
    }

    // Equilibrium concentration at T
    const Real EmpiricalEquilibrium(Real T)
    {
        return pow(10.0, analytical_expression[0] + analytical_expression[1] * T + analytical_expression[2] / T + analytical_expression[3] * log10(T));
    }

    ///A
    const Real DebyeHuckelParam(Real T)
    {
        Real A = sqrt(2 * pi * N_A * rho_w) * pow((e * e) / (4 * pi * eps_0 * eps_r * k_b * T), 1.5); //kg/mol
        return A;
    }

    const Real DebyeHuckelActivity(Real T, Real yEtc1, Real yEtc2)
    {
        return exp(-DebyeHuckelParam(T) * (Z_A * Z_A) * IonicStrength(yEtc1, yEtc2));
    }

    const Real ExtendedDebyeHuckelActivity(Real T, Real yEtc1, Real yEtc2)
    {
        Real I = IonicStrength(yEtc1, yEtc2);
        return exp(-DebyeHuckelParam(T) * fabs(Z_A * Z_B) * sqrt(I) / (1 + ExtendedDebyeHuckelParam(T) * sqrt(I)));
    }

    const Real ExtendedDebyeHuckelParam(Real T)
    {
        return sqrt(2 * e * e * N_A * rho_w / (eps_0 * eps_r * k_b * T));
    }

    // Activity coefficient from Pitzer's eq.
    const Real pitzerActivityCoefficient(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return pitzerActivityCoefficient(T, IonicStrength(yEtc1, yEtc2), MeanMolality(yA, yB, yEtc1, yEtc2));
    }

    // Activity coefficient from Pitzer's eq.
    const Real pitzerActivityCoefficient(Real T, Real I, Real meanMolality)
    {
        if (I < SMALL)
            return 1;

        const Real alpha_1 = 1.4;
        const Real alpha_2 = 12;
        const Real b = 1.2;

        const Real A = DebyeHuckelParam(T); //kg/mol
        const Real B_gamma = 2 * beta_0 + 2 * beta_1 / ((alpha_1 * alpha_1) * I) * (1 - (1 + alpha_1 * sqrt(I) - 0.5 * (alpha_1 * alpha_1) * I) * exp(-alpha_1 * sqrt(I))) + 2 * beta_2 / ((alpha_2 * alpha_2) * I) * (1 - (1 + alpha_2 * sqrt(I) - 0.5 * (alpha_2 * alpha_2) * I) * exp(-alpha_2 * sqrt(I)));

        const Real f_gamma = -A / 3 * (sqrt(I) / (1 + b * sqrt(I)) + 2 / b * log(1 + b * sqrt(I)));

        const Real C_gamma = 3 / 2 * C_Phi;

        const Real ln_gamma = fabs(Z_A * Z_B) * f_gamma + meanMolality * (2 * nu_A * nu_B / nu) * B_gamma + pow(meanMolality, 2) * (2 * pow(nu_A * nu_B, 1.5) / nu) * C_gamma;

        return exp(ln_gamma);
    }

    // Nucleation rate
    const Real NucleationRate(Real S, Real T)
    {
        if (S <= 1.001)
        {
            return 0;
        }

        const Real Vm = MolarMass / rho;
        const Real ExponentNumerator = 16.0 * N_A * pi * pow(Vm, 2.0) * pow(sigma, 3.0);
        const Real ExponentDenominator = 3.0 * pow(R * T, 3.0) * pow(log(S), 2.0);

        return An * exp(-ExponentNumerator / ExponentDenominator);
    }

    // Nucleation Critical Radius
    const Real CriticalRadius(Real S, Real T)
    {
        if (S <= 1.001)
        {
            return 0;
        }

        const Real ExponentNumerator = 2.0 * MolarMass * sigma;
        const Real ExponentDenominator = rho * R * T * log(S);

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

    // Saturation Equation
    const Real ExtendedDebyeHuckelSaturationRate(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return pow(MeanMolality(yA, yB, yEtc1, yEtc2) * ExtendedDebyeHuckelActivity(T, yEtc1, yEtc2), nu) / Equilibrium(T);
    }

    // Saturation Equation
    const Real ExtendedDebyeHuckelSaturationIndex(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return log10(ExtendedDebyeHuckelSaturationRate(T, yA, yB, yEtc1, yEtc2));
    }

    // Saturation Equation
    const Real DebyeHuckelSaturationRate(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return pow(MeanMolality(yA, yB, yEtc1, yEtc2) * DebyeHuckelActivity(T, yEtc1, yEtc2), nu) / Equilibrium(T);
    }

    // Saturation Equation
    const Real DebyeHuckelSaturationIndex(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return log10(DebyeHuckelSaturationRate(T, yA, yB, yEtc1, yEtc2));
    }

    // Iterative saturation to avoid overreaction - not in use
    const Real IterativePitzerSaturationIndex(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return log10(IterativePitzerSaturationRate(T, yA, yB, yEtc1, yEtc2));
    }

    // Wall saturation model
    const Real PitzerWallConcentrationA(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        double SR = IterativePitzerSaturationRate(T, yA, yB, yEtc1, yEtc2);
        double mA = TotalMolality(yEtc1,yEtc2)*yA;
        return fmin(mA, 1/SR * mA);
    }

    // Wall saturation model
    const Real PitzerWallConcentrationB(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        double SR = IterativePitzerSaturationRate(T, yA, yB, yEtc1, yEtc2);
        double mB = TotalMolality(yEtc1,yEtc2)*yB;
        return fmin(mB, 1/SR * mB);
    }

    // Compute Ionic Strength
    const Real IonicStrength(Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        const Real m[2] = {yEtc1 * mTot, yEtc2 * mTot};
        const Real Z[2] = {1, 2};
        return ChemistryFunctions::IonicStrength(m, Z, 2);
    }

    // Solve for saturation correcting for change in activity on reaction
    const Real IterativePitzerSaturationRate(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        Real mA = yA * mTot;
        Real mB = yB * mTot;

        if (fmin(mA, mB) < SMALL)
            return SMALL;

        const int noCorrectionSteps = 6;
        const Real relaxationFactor = 0.8; // Important due to mathematical incorrectness of correction step

        const Real equilibrium = Equilibrium(T);
        const Real meanMolality = MeanMolality(mA, mB);
        const Real I = IonicStrength(yEtc1, yEtc2);
        const Real gamma = pitzerActivityCoefficient(T, I, meanMolality);

        Real S;        // Resulting saturation index
        Real Sc = 0;   // Corrected
        Real mSat = 0; // Oversaturated molality
        
        S = pow(meanMolality * gamma, nu) / equilibrium;

        return S;
        
        for (int i = 0; i < noCorrectionSteps + 1; i++)
        {

            // Correcting for change in ionic strength
            Real newI = fmax(SMALL, I - 0.5 * mSat * (Z_A * Z_A + Z_B * Z_B));
            Real newGamma = pitzerActivityCoefficient(T, newI, meanMolality);

            //return newGamma;
            Real m_c = fmax(SMALL, meanMolality - mSat); // Corrected m -> K
            // Sc = SaturationRate(_equilibrium, newGamma, m_c); // Corrected Saturation --->1
            Sc = pow(m_c * newGamma, nu) / equilibrium;

            if (Sc > 1)
            {
                // mSat = mSat + relaxationFactor * (Sc - 1 + SMALL)/(SC + SMALL) * (fmin(mA, mB) - mSat);
            }
            else
            {
                mSat = 0;
            }

            // S = SaturationRate(_equilibrium, newGamma, _MeanMolality);
            S = pow(meanMolality * newGamma, nu) / equilibrium;
        }
        return S;
    }

private:
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
            
        return mMean * nu_i / (nu * m_i);
    }

    // SaturationRate to activity
    const Real dSR_dgamma(Real mMean, Real K)
    {
        // TODO: find derivative
        return mMean / K;
    }

    // activity to ionicStrength
    const Real dgamma_dI_debyeHuckel(Real Temperature, Real Z_i)
    {
        return Z_i * Z_i * ExtendedDebyeHuckelParam(Temperature);
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

    // activity to mole fraction
    const Real dgamma_dyi_DebyeHuckel(Real Temperature, Real y_i, Real Z_i, Real mTot)
    {
        return dgamma_dI_debyeHuckel(Temperature, Z_i) * dI_dyi(y_i, Z_i, mTot);
    }

    const Real dS_dyi_DebyeHuckel(Real Temperature, Real y_i, Real m_i, Real nu_i, Real Z_i, Real gamma, Real K, Real mMean, Real mTot)
    {
        return dSR_dmMean(gamma, K) * dmMean_dmi(mMean, m_i, nu_i) * dm_dyi(y_i, mTot) + dSR_dgamma(mMean, K) * dgamma_dyi_DebyeHuckel(Temperature, y_i, Z_i, mTot); // not done
    }

    const Real dRw_dyi(Real Temperature, Real y_i, Real m_i, Real nu_i, Real Z_i, Real mMean, Real SaturationRate, Real gamma, Real K, Real mTot, Real drivingMolality)
    {
        // TODO derivatibe of drivingMolalityRate?
        return DrivingMolalityRate(m_i, drivingMolality) * dm_dyi(y_i, mTot) + dRw_dSR(SaturationRate, drivingMolality) * dS_dyi_DebyeHuckel(Temperature, y_i, m_i, nu_i, Z_i, gamma, K, mMean, mTot);
    }
};

#endif // CHEMISTRY_H
