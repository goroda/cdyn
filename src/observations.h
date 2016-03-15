// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <stdlib.h>

struct Observation;
struct Observation * 
observation_alloc(size_t, size_t,
                  int (*)(double,const double*,double*,void*),
                  void *);
void observation_free(struct Observation *);
size_t observation_get_dx(const struct Observation *);
size_t observation_get_dy(const struct Observation *);
int observation_observe(const struct Observation *,
                        double, const double *, double *);
#endif
