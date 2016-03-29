// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "dynamics.h"


struct Integrator;
struct Integrator *
integrator_create(size_t,
                  int (*)(double,const double*,double*,double*,void*),
                  void *);
void integrator_destroy(struct Integrator *);
void integrator_set_type(struct Integrator *, char *);
void integrator_set_dt(struct Integrator *, double);
void integrator_set_adaptive_opts(struct Integrator *, double, double, double);
void integrator_set_verbose(struct Integrator *, int);
void integrator_step(struct Integrator *,
                     double, double,
                     const double *, double *);


int euler_step(double time, const double * x, double * nextx,
                   double dt, struct Dyn * dyn, double * driftv);

int rk4_step(double, const double *, double *,
             double,struct Dyn *, double *);
int 
euler_maruyama_step(double, const double *, const double *,
                        double *, double, struct Dyn *, 
                        double *, double *);

#endif
