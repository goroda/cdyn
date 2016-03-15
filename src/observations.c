// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#include <stdlib.h>
#include <stdio.h>
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
    return h;
}

/*!
  Free observation operator
*/
void observation_free(struct Observation * obs)
{
    if (obs != NULL){
        free(obs); obs = NULL;
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
  Generate an observation
*/
int observation_observe(const struct Observation * obs,
                        double time, const double * x, double * y)
{

    assert (obs != NULL);
    int res = obs->obs(time,x,y,obs->obsarg);
    return res;
}
