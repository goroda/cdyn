// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <assert.h>

#include "observations.h"

/** \struct Observation
 *  \brief Observation operator
 *  \var Observation::dx
 *  dimension of state space
 *  \var Observation::dy
 *  dimension of observation
 *  \var Observation::obs
 *  observation dynamics
 *  f(time,state,out,args)
 *  \var Observation::obsarg
 *  Additional arguments to observations
 */
struct Observation
{
    size_t dy;
    size_t dx;

    int (*obs)(double,const double *,double *,void *);

    double * cov;
    void * obsarg;
};

/*!
  Allocate observation operator

  \param[in] dx   - dimension of state
  \param[in] dy   - dimension of observation
  \param[in] obs  - observation function
  \param[in] oarg - additional arguments to observation function

  \return observation operator
*/
struct Observation * 
observation_alloc(size_t dx, size_t dy,
                  int (*obs)(double,const double*,double*,void*),
                  void * oarg)
{
    struct Observation * h = malloc(sizeof(struct Observation));
    if (h == NULL){
        fprintf(stderr,"Error allocating memory for Observation\n");
        exit(1);
    }

    h->dx = dx;
    h->dy = dy;
    
    h->obs = obs;
    h->obsarg = oarg;

    h->cov = NULL;
    return h;
}

/*!
  Free observation operator
*/
void observation_free(struct Observation * obs)
{
    if (obs != NULL){
        free(obs->cov); obs->cov = NULL;
        free(obs); obs = NULL;
    }
}

/*!
  Add noise covariance
*/
void observation_set_noise_cov(struct Observation * obs, const double * cov)
{
    if (obs != NULL){
        free(obs->cov);
        obs->cov = malloc(obs->dy*obs->dy * sizeof(double));
        if (obs->cov == NULL){
            fprintf(stderr,"Cannot allocate cov in observation_set_noise_cov\n");
            exit(1);
        }
        memmove(obs->cov,cov, obs->dy*obs->dy * sizeof(double));
    }
}

/*!
  Get the dimension of state space
*/
size_t observation_get_dx(const struct Observation * obs)
{
    assert (obs != NULL);
    return obs->dx;
}

/*!
  Get the dimension of the observations
*/
size_t observation_get_dy(const struct Observation * obs)
{
    assert (obs != NULL);
    return obs->dy;
}

/*!
  Get the noise covariance
*/
double * observation_get_noise_cov(const struct Observation * obs)
{
    assert (obs != NULL);
    return obs->cov;
}

/*!
  Generate an observation
*/
int observation_observe(const struct Observation * obs,
                        const double * noise,
                        double time, const double * x, double * y)
{

    assert (obs != NULL);
    int res = obs->obs(time,x,y,obs->obsarg);
    if (noise != NULL){
        for (size_t ii = 0; ii < obs->dy; ii++){
            y[ii] += noise[ii];
        }
    }
    return res;
}

///////////////////////////////////////////////////////
// Interface stuff for FT

/*!
  Initialize  observation and time couple
*/
void observation_plus_time_init(struct ObservationPlusTime * ot,
                                const struct Observation * obs,
                                double time)
{
    ot->obs = obs;
    ot->time = time;
}

/*!
  Set time of observation and time couple
*/
void observation_plus_time_set_time(struct ObservationPlusTime * ot,
                                    double time)
{
    ot->time = time;
}

/*!
  Interface for approximation by ft
*/
double observation_ft_bb(double * x, size_t ind, void * arg)
{

    struct ObservationPlusTime * ot = arg;
    size_t dy = ot->obs->dy;
    double * out = malloc(dy*sizeof(double));
    if (out == NULL){
        fprintf(stderr,"Memory error in observation_ft_bb\n");
        exit(1);
    }
    int res = ot->obs->obs(ot->time,x,out,ot->obs->obsarg);
    assert (res == 0);
    double eval = out[ind];
    free(out); out = NULL;
    return eval;
}
