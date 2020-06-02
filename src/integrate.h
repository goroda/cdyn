// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <stdlib.h>

struct Integrator;
struct Integrator *
integrator_create(size_t,
                  int (*)(double,const double*,double*,double*,void*),
                  void *);
struct Integrator *
integrator_create_controlled(
    size_t, size_t,
    int (*)(double, const double *, const double *,double*,double *, void *),
    void *,
    int (*)(double, const double *, double *, void *),
    void *);

int integrator_update_control_args(struct Integrator * i, void * args);

int integrator_eval_controller(struct Integrator *,
                               double, const double *, double *);
void integrator_destroy(struct Integrator *);
void integrator_set_type(struct Integrator *, char *);
void integrator_set_dt(struct Integrator *, double);
void integrator_set_adaptive_opts(struct Integrator *, double, double, double);
void integrator_set_verbose(struct Integrator *, int);
void integrator_set_dargs(struct Integrator *, void *);
void integrator_step(struct Integrator *,
                     double, double,
                     const double *, double *);

/* int  */
/* euler_maruyama_step(double, const double *, const double *, */
/*                         double *, double, struct Dyn *,  */
/*                         double *, double *); */

#endif
