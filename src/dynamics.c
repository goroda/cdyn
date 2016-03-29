// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "dynamics.h"

/** \struct Drift
 *  \brief Drift dynamics
 *  \var Drift::dx
 *  dimension of state space
 *  \var Drift::b
 *  uncontrolled dynamics RHS of drift term to stochastic differential equation
 *  f(time,state,jac,out,args)
 *  \var Drift::bargs
 *  Additional arguments to dynamics
 */
struct Drift
{

    size_t dx;
    int (*b)(double,const double *,double *,double *,void *);
    void * bargs;
};

/*!
  Add uncontrolled drift dynamcis

  \param[in,out] dr    - allocated drift object
  \param[in]     b     - uncontrolled dynamics
  \param[in]     bargs - arguments to dynamics

*/
void drift_add_rhs(struct Drift * dr,
                   int (*b)(double,const double*,double*,double*,void*),
                   void * bargs)
{
    assert (dr != NULL);
    dr->b = b;
    dr->bargs = bargs;
}

/*!
  Allocate uncontrolled drift dynamics

  \param[in] dx - dimension of state
  \param[in] bc - uncontrolled dynamics
  \param[in] ba - arguments to dynamics

  \return drift dynamics
*/
struct Drift * 
drift_alloc(size_t dx,
            int (*bc)(double,const double*,double *,double*,void*),
            void * ba)
{
    struct Drift * b = malloc(sizeof(struct Drift));
    if (b == NULL){
        fprintf(stderr,"Error allocating memory for Drift\n");
        exit(1);
    }

    b->dx = dx;
    drift_add_rhs(b,bc,ba);
    return b;
}


/*!
  Free drift dynamics
*/
void drift_free(struct Drift * drift)
{
    if (drift != NULL){
        free(drift); drift = NULL;
    }
}

/*!
  Get state space size
*/
size_t drift_get_dx(const struct Drift * b)
{
    assert (b != NULL);
    return b->dx;
}

/*!
  Evaluate the dynamics
*/
int drift_eval(const struct Drift * b, double time, 
               const double * x, double * out,double * jac)
{
    int res;
    assert ( b != NULL);
    res = b->b(time,x,out,jac,b->bargs);
    return res;
}

////////////////////////////////////////////////////////////////////////

/** \struct Diff
 *  \brief Diffusion dynamics
 *  \var Diff::dx
 *  dimension of state space
 *  \var Diff::dw
 *  dimension of random walk
 *  \var Diff::s
 *  RHS of diffusion term of stochastic differential equation
 *  f(time,state,out,args)
 *  \var Diff::sargs
 *  Additional arguments to dynamics
 */
struct Diff
{
    size_t dx;
    size_t dw;

    int (*s)(double,const double*,double*,void*);
    void * sargs;
};

/*!
  Add uncontrolled drift dynamcis

  \param[in,out] diff  - allocated diff object
  \param[in]     s     - uncontrolled dynamics
  \param[in]     sargs - arguments to dynamics

*/
void diff_add_rhs(struct Diff * diff,
                  int (*s)(double,const double*,double*,
                           void*),
                  void * sargs)
{
    assert (diff != NULL);
    diff->s = s;
    diff->sargs = sargs;
}

/*!
  Allocate uncontrolled diffusion dynamics

  \param[in] dx     - dimension of state
  \param[in] dw     - dimension of noise
  \param[in] sd     - uncontrolled dynamics
  \param[in] sdargs - arguments to dynamics

  \return diffusion dynamics
*/
struct Diff * 
diff_alloc(size_t dx, size_t dw,
           int (*sd)(double,const double*,double*,void*),
           void * sdargs)
{
    struct Diff * s = malloc(sizeof(struct Diff));
    if (s == NULL){
        fprintf(stderr,"Error allocating memory for Diff\n");
        exit(1);
    }
    
    s->dx = dx;
    s->dw = dw;
    diff_add_rhs(s,sd,sdargs);

    return s;
}

/*!
  Free diffusion dynamics
*/
void diff_free(struct Diff * diff)
{
    if (diff != NULL){
        free(diff); diff = NULL;
    }
}

/*!
   Evaluate diffusion dynamics
*/
int diff_eval(const struct Diff * b, double time, const double * x, 
              double * out)
{
    int res;
    assert(b->s != NULL);
    res = b->s(time,x,out,b->sargs);
    
    return res;
}

/*!
   Get the dimension of the noise
*/
size_t diff_get_dw(struct Diff * b){
    assert (b != NULL);
    return b->dw;
}

/////////////////////////////////////////////////////////

/** \struct Dyn
 *  \brief Stochastic Differential Equation Dynamics
 *  \var Dyn::drift
 *  drift dynamics
 *  \var Dyn::diff
 *  diffusion dynamics
 *  \var Dyn::s
 *  RHS of diffusion term of stochastic differential equation
 *  \var Dyn::sargs
 *  Additional arguments to dynamics
 */
struct Dyn
{
    struct Drift * drift;
    struct Diff  * diff;
};


/*!
   Create dynamic object
*/
struct Dyn * dyn_alloc(struct Drift * drift, struct Diff * diff)
{
    struct Dyn * dyn = malloc(sizeof(struct Dyn));
    if (dyn == NULL){
        fprintf(stderr,"Cannot allocate memory for Dynamics\n");
        exit(1);
    }
    dyn->drift = drift;
    dyn->diff = diff;

    return dyn;
}

/*!
   Free dynamics object
*/
void dyn_free(struct Dyn * dyn)
{
    if (dyn != NULL){
        free(dyn); dyn = NULL;
    }
}

/*!
   Free dynamics object along with drift and diffusion dynamics
*/
void dyn_free_deep(struct Dyn * dyn)
{
    if (dyn != NULL){
        drift_free(dyn->drift); dyn->drift = NULL;
        diff_free(dyn->diff); dyn->diff = NULL;
        dyn_free(dyn); dyn = NULL;
    }
}

/*!
   Get size of state space
*/
size_t dyn_get_dx(const struct Dyn * dyn)
{
    assert (dyn != NULL);
    return drift_get_dx(dyn->drift);
}

/*!
   Get size of stochastic space
*/
size_t dyn_get_dw(const struct Dyn * dyn)
{
    assert (dyn != NULL);
    return diff_get_dw(dyn->diff);
}

/*!
   Evaluate dynamics
*/
int dyn_eval(const struct Dyn * dyn,double time,
             const double * x,double * drift, 
             double * jacdr,double * diff)
{
    int res = 0;
    if (drift != NULL){
        res = drift_eval(dyn->drift,time,x,drift,jacdr);
    }
    if (res != 0){
        return res;
    }
    if (diff != NULL){
        res = diff_eval(dyn->diff,time,x,diff);
    }
    return res;
}


