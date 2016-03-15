// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "dynamics.h"
int euler_step(double time, const double * x, double * nextx,
                   double dt, struct Dyn * dyn, double * driftv);
int 
euler_maruyama_step(double, const double *, const double *,
                        double *, double, struct Dyn *, 
                        double *, double *);

#endif
