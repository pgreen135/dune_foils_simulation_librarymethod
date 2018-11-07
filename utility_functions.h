#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H

#include <vector>
#include "TVector3.h"

//A large number of these are simple functions needed to create the distributions
//such as beta decay or a poisson distribution.
//The other half of the code deals with parameterisations that Diego created
//to express the time it takes photons to propagate through the detector.

namespace utility{

    int poisson(double mean, double draw, double eng);
    double SpectrumFunction(double *x, double *par);
    double fsn(double *x, double *par);
    double fso(double *x, double *par);
    double Rn_function(double *x, double *par);
    double Scintillation_function(double *t, double *par);
    
  }

#endif
