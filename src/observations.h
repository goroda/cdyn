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
void observation_set_noise_cov(struct Observation *, const double *);
size_t observation_get_dx(const struct Observation *);
size_t observation_get_dy(const struct Observation *);
double * observation_get_noise_cov(const struct Observation *);
int observation_observe(const struct Observation *, const double *, 
                        double, const double *, double *);

struct ObservationPlusTime
{
    const struct Observation * obs;
    double time;
};

void observation_plus_time_init(struct ObservationPlusTime *,
                                const struct Observation *,
                                double);
void observation_plus_time_set_time(struct ObservationPlusTime *,
                                    double);
double observation_ft_bb(double *, size_t, void *);

#endif
