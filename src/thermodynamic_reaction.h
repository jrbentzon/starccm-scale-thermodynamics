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

class ThermodynamicReaction
{
public:
    SimpleReaction &reaction;
    EquilibriumFormulation &equilibriumFormulation;
    ActivityModel &activityModel;

    // numerical
    const Real SMALL = 1e-16;

    // Constructor
    ThermodynamicReaction(SimpleReaction &reaction,
                          EquilibriumFormulation &equilibriumFormulation,
                          ActivityModel &activityModel)
        : reaction(reaction),
          activityModel(activityModel),
          equilibriumFormulation(equilibriumFormulation)
    {
    }

    // Mean Molarity
    const void MeanMolality(Real *meanMolality, int size, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            meanMolality[i] = reaction.MeanMolality(yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }

    // Equilibrium Function
    const void Equilibrium(Real *equilibrium, int size, Real *Temperature)
    {
        for (int i = 0; i < size; i++)
        {
            equilibrium[i] = equilibriumFormulation.Equilibrium(Temperature[i]);
        }
    }

    // Pitzers activitity coeff
    const void ActivityCoefficient(Real *gamma, int size, Real *Temperature, Real *yA, Real *yB, Real *yEtc1, Real *yEtc2)
    {
        for (int i = 0; i < size; i++)
        {
            gamma[i] = activityModel.ActicityCoefficient(Temperature[i], yA[i], yB[i], yEtc1[i], yEtc2[i]);
        }
    }
};

#endif // BARITEREACTION_H
