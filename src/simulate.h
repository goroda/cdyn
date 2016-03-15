// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#ifndef SIMULATE_H
#define SIMULATE_H

#include <stdlib.h>
#include <stdio.h>

#include "dynamics.h"

struct Trajectory;
struct Trajectory * trajectory_init(size_t, size_t, double, 
                                    const double *, 
                                    const double *);
void trajectory_free(struct Trajectory *);
int trajectory_add(struct Trajectory **, size_t, size_t,
                   double,const double *,const double *);
void trajectory_print(const struct Trajectory *, FILE *, int);

double trajectory_get_last_time(const struct Trajectory *);
double * trajectory_get_last_state(const struct Trajectory *);
double * trajectory_get_last_control(const struct Trajectory *);

int trajectory_step(struct Trajectory *,
                    struct Dyn *,
                    double, char *,
                    double *,
                    void *, void *);

#endif  
