#ifndef _CHI2_TEST_HPP_
#define _CHI2_TEST_HPP_
#include <string>
#include "ISRSolverSLAE.hpp"

void chi2TestModel(int, double,
                   ISRSolverSLAE*,
                   const std::string&,
                   const std::string&,
                   const std::string&,
                   const std::string&);

void chi2TestData(int, double,
                  ISRSolverSLAE*,
                  const std::string&);

#endif
