#ifndef _CHI2_TEST_HPP_
#define _CHI2_TEST_HPP_
#include <string>
#include "ISRSolverSLE.hpp"

// !!! TO DO : keep it simple ...
void chi2TestModel(int, double,
                   ISRSolverSLE*,
                   const std::string&,
                   const std::string&,
                   const std::string&,
                   const std::string&);

void chi2TestData(int, double,
                  ISRSolverSLE*,
                  const std::string&);

#endif
