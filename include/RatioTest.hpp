#pragma once

#include "ISRSolverStructs.hpp"
#include "ISRSolverSLE.hpp"

/**
 * This function is used to get ratio of a numerical solution to
 a model Born cross section averaged over a number of numerical experiments
 * @param solver a solver
 * @param args an input arguments
 * @see ISRSolverSLE
 * @see RatioTestModelArgs
 */
void ratioTestModel(ISRSolverSLE* solver,
                    const RatioTestModelArgs& args);
