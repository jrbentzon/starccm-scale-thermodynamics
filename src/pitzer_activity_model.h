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

#ifndef PITZER_ACTIVITY_H
#define PITZER_ACTIVITY_H

//#include "chemistry.h"
#include "uclib.h"
#include "math.h"
#include <cstdlib>
#include "activity_model.h"
#include "simple_reaction.h"

class PitzerActivityModel : public ActivityModel
{
public:
    // Pitzer model parameters (input)
    Real beta_0;
    Real beta_1;
    Real beta_2;
    Real C_Phi;

    // Reaction
    SimpleReaction &reaction;

public:
    PitzerActivityModel(SimpleReaction &reaction,
                        Real beta_0,
                        Real beta_1,
                        Real beta_2,
                        Real C_Phi)
        : beta_0(beta_0),
          beta_1(beta_1),
          beta_2(beta_2),
          C_Phi(C_Phi),
          reaction(reaction)
    {
    }

    Real getActivityCoefficient(Real T, Real yA, Real yB, Real yEtc1, Real yEtc2)
    {
        return IonicStrength(yEtc1, yEtc2);
        //return pitzerActivityCoefficient(T, IonicStrength(yEtc1, yEtc2), reaction.MeanMolality(yA, yA, yEtc1, yEtc2));
    };

    // Activity coefficient from Pitzer's eq.
    Real pitzerActivityCoefficient(Real T, Real I, Real meanMolality)
    {
        // Pitzer model constants
        Real alpha_1 = 1.4;
        Real alpha_2 = 12;
        Real b = 1.2;

        Real SMALL = 1e-16;

        if (I < SMALL)
            return 1.0;

        Real A = DebyeHuckelParam(T); //kg/mol
        Real B_gamma = 2 * beta_0 + 2 * beta_1 / ((alpha_1 * alpha_1) * I) * (1 - (1 + alpha_1 * sqrt(I) - 0.5 * (alpha_1 * alpha_1) * I) * exp(-alpha_1 * sqrt(I))) + 2 * beta_2 / ((alpha_2 * alpha_2) * I) * (1 - (1 + alpha_2 * sqrt(I) - 0.5 * (alpha_2 * alpha_2) * I) * exp(-alpha_2 * sqrt(I)));

        Real f_gamma = -A / 3 * (sqrt(I) / (1 + b * sqrt(I)) + 2 / b * log(1 + b * sqrt(I)));

        Real C_gamma = 3 / 2 * C_Phi;

        Real ln_gamma = fabs(reaction.Z_A * reaction.Z_B) * f_gamma + meanMolality * (2 * reaction.nu_A * reaction.nu_B / reaction.nu()) * B_gamma + pow(meanMolality, 2) * (2 * pow(reaction.nu_A * reaction.nu_B, 1.5) / reaction.nu()) * C_gamma;

        return exp(ln_gamma);
    };

    ///A
    Real DebyeHuckelParam(Real T)
    {
        Real pi = M_PI;
        Real rho_w = 997;
        Real N_A = 6.022140857e23;
        Real e = 1.60206e-19;
        Real eps_0 = 8.8542e-12;
        Real eps_r = 78.4;
        Real k_b = 1.38064852e-23;

        Real A = sqrt(2 * pi * N_A * rho_w) * pow((e * e) / (4 * pi * eps_0 * eps_r * k_b * T), 1.5); //kg/mol
        return A;
    };

    // Compute Ionic Strength
    Real IonicStrength(Real yEtc1, Real yEtc2)
    {
        return 1;
        Real mTot = reaction.TotalMolality(yEtc1, yEtc2);
        Real m[2] = {yEtc1 * mTot, yEtc2 * mTot};
        Real Z[2] = {1, 2};
        Real I = 0;
        for (int j = 0; j < 2; j++)
        {
            I = I + m[j] * pow(Z[j], 2); //0.5 should be taken outside loop for better performance
        }
        I = 0.5 * I;
        return I;
    };
};

#endif // PITZER_ACTIVITY_H
