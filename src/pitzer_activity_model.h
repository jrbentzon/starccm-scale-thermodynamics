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

#ifndef PITZER_ACTIVITY_MODEL_H
#define PITZER_ACTIVITY_MODEL_H

#include "chemistry.h"
#include "uclib.h"
#include "math.h"
#include <cstdlib>
#include "activity_model.h"
#include "simple_reaction.h"

class PitzerActivityModel : public ActivityModel
{
public:
    const Real &SMALL = 1e-16;
    const Real &N_A = ChemistryFunctions::N_A();
    const Real &rho_w = ChemistryFunctions::densityWater();
    const Real &e = ChemistryFunctions::electronicCharge();
    const Real &eps_0 = ChemistryFunctions::permittivityVacuum();
    const Real &eps_r = ChemistryFunctions::permittivityWater();
    const Real &k_b = ChemistryFunctions::k_b();

    // Pitzer model parameters (input)
    const Real &beta_0;
    const Real &beta_1;
    const Real &beta_2;
    const Real &C_Phi;

    // Reaction
    const Real &nu_A;
    const Real &nu_B;
    const Real &Z_A;
    const Real &Z_B;

    PitzerActivityModel(const Real &nu_A,
                        const Real &nu_B,
                        const Real &Z_A,
                        const Real &Z_B,
                        const Real &beta_0,
                        const Real &beta_1,
                        const Real &beta_2,
                        const Real &C_Phi) : beta_0(beta_0),
                                             beta_1(beta_1),
                                             beta_2(beta_2),
                                             C_Phi(C_Phi),
                                             nu_A(nu_A),
                                             nu_B(nu_B),
                                             Z_A(Z_A),
                                             Z_B(Z_B)
    {
    }

    // Activity Coeffiecient (gamma)
    const Real ActivityCoefficient(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return pitzerActivityCoefficient(T, IonicStrength(yEtc1, yEtc2), MeanMolality(yA, yB, yEtc1, yEtc2));
    }

    // Activity coefficient from Pitzer's eq.
    const Real pitzerActivityCoefficient(Real T, Real I, Real meanMolality)
    {
        const Real alpha_1 = 1.4;
        const Real alpha_2 = 12;
        const Real b = 1.2;

        if (I < SMALL)
            return 1;

        const Real A = DebyeHuckelParam(T); //kg/mol
        const Real B_gamma = 2 * beta_0 + 2 * beta_1 / ((alpha_1 * alpha_1) * I) * (1 - (1 + alpha_1 * sqrt(I) - 0.5 * (alpha_1 * alpha_1) * I) * exp(-alpha_1 * sqrt(I))) + 2 * beta_2 / ((alpha_2 * alpha_2) * I) * (1 - (1 + alpha_2 * sqrt(I) - 0.5 * (alpha_2 * alpha_2) * I) * exp(-alpha_2 * sqrt(I)));

        const Real f_gamma = -A / 3 * (sqrt(I) / (1 + b * sqrt(I)) + 2 / b * log(1 + b * sqrt(I)));

        const Real C_gamma = 3 / 2 * C_Phi;

        const Real ln_gamma = fabs(Z_A * Z_B) * f_gamma + meanMolality * (2 * nu_A * nu_B / (nu_A + nu_B)) * B_gamma + pow(meanMolality, 2) * (2 * pow(nu_A * nu_B, 1.5) / (nu_A + nu_B)) * C_gamma;

        return exp(ln_gamma);
    }

    ///A
    const Real DebyeHuckelParam(Real T)
    {
        const Real &rho_w = 997;
        return sqrt(2 * M_PI * 6.022140857e23 * rho_w) * pow(2.131077129519e-7 / T, 1.5); //kg/mol
    }

    // Compute Ionic Strength
    const Real IonicStrength(Real yEtc1, Real yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        const Real m[2] = {yEtc1 * mTot, yEtc2 * mTot};
        const Real Z[2] = {1, 2};
        return ChemistryFunctions::IonicStrength(m, Z, 2);
    }

    // Total Concentration
    const Real TotalMolality(Real &yEtc1, Real &yEtc2)
    {
        return pow(ChemistryFunctions::MolarMassOfWater() / (1 - (yEtc1 + yEtc2)) + SMALL, -1);
    }

    // Returns Mean molality
    const Real MeanMolality(Real &yA, Real &yB, Real &yEtc1, Real &yEtc2)
    {
        const Real mTot = TotalMolality(yEtc1, yEtc2);
        return MeanMolality(yA * mTot, yB * mTot);
    }

    // Mean Concentration
    const Real MeanMolality(Real mA, Real mB)
    {
        return pow(pow(fmax(mA, SMALL), nu_A) * pow(fmax(mB, SMALL), nu_B) + SMALL, 1.0 / (nu_A + nu_B));
    }
};

#endif // PITZER_ACTIVITY_MODEL_H
