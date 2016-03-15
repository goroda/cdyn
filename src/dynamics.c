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
 *  \var Drift::controlled
 *  1 for controlled system 0 for no
 *  \var Drift::dx
 *  dimension of state space
 *  \var Drift::du
 *  dimension of control
 *  \var Drift::bc
 *  controlled dynamics RHS of drift term to stochastic differential equation
 *  f(time,state,control,out,grad,args) (gradient of control)
 *  \var Drift::bargs
 *  Additional arguments to dynamics
 *  \var Drift::b
 *  uncontrolled dynamics RHS of drift term to stochastic differential equation
 *  f(time,state,out,args)
 *  \var Drift::bargs
 *  Additional arguments to dynamics
 */
struct Drift
{

    int controlled;

    size_t dx;
    size_t du;

    int (*bc)(double,const double *,const double *,double *,double*,void *);
    void * bcargs;

    int (*b)(double,const double *,double *,void *);
    void * bargs;
};

/*!
  Add controlled drift dynamcis

  \param[in,out] dr    - allocated drift object
  \param[in]     b     - controlled dynamics
  \param[in]     bargs - arguments to dynamics

*/
void drift_add_rhs_controlled(struct Drift * dr,
                              int (*b)(double,const double*,const double*,
                                       double*,double*,void*),
                              void * bargs)
{
    assert (dr != NULL);
    dr->bc = b;
    dr->bcargs = bargs;
}

/*!
  Allocate controlled drift dynamics

  \param[in] dx - dimension of state
  \param[in] du - dimension of control
  \param[in] bc - controlled dynamics
  \param[in] ba - arguments to dynamics

  \return drift dynamics
*/
struct Drift * drift_alloc_controlled(size_t dx, size_t du,
                                      int (*bc)(double,const double*,const double*,
                                               double*,double*,void*),
                                      void * ba)
{
    struct Drift * b = malloc(sizeof(struct Drift));
    if (b == NULL){
        fprintf(stderr,"Error allocating memory for Drift\n");
        exit(1);
    }

    b->controlled = 1;
    b->dx = dx;
    b->du = du;
    
    drift_add_rhs_controlled(b,bc,ba);
    b->b = NULL;
    b->bargs = NULL;

    return b;
}

/*!
  Add uncontrolled drift dynamcis

  \param[in,out] dr    - allocated drift object
  \param[in]     b     - uncontrolled dynamics
  \param[in]     bargs - arguments to dynamics

*/
void drift_add_rhs(struct Drift * dr,
                   int (*b)(double,const double*,double*,void*),
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
drift_alloc_uncontrolled(size_t dx,
                         int (*bc)(double,const double*,double*,void*),
                         void * ba)
{
    struct Drift * b = malloc(sizeof(struct Drift));
    if (b == NULL){
        fprintf(stderr,"Error allocating memory for Drift\n");
        exit(1);
    }

    b->controlled = 0;
    b->dx = dx;
    b->du = 0;

    b->bc = NULL;
    b->bcargs = NULL;
    
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
  Get control space size
*/
size_t drift_get_du(const struct Drift * b)
{
    assert (b != NULL);
    return b->du;
}

/*!
  Evaluate the dynamics
*/
int drift_eval(const struct Drift * b, double time, 
               const double * x,
               const double * u, double * out, double * jac)
{
    int res;
    assert ( b != NULL);
    if (b->controlled == 0){
        assert(b->b != NULL);
        /* assert(u == NULL); */
        /* assert(jac == NULL); */
        res = b->b(time,x,out,b->bargs);
    }
    else{
        assert(b->bc != NULL);
        res = b->bc(time,x,u,out,jac,b->bcargs);
    }

    return res;
}

////////////////////////////////////////////////////////////////////////

/** \struct Diff
 *  \brief Diffusion dynamics
 *  \var Diff::dx
 *  dimension of state space
 *  \var Diff::du
 *  dimension of control
 *  \var Diff::dw
 *  dimension of random walk
 *  \var Diff::sc
 *  RHS of controlled diffusion term of stochastic differential equation
 *  f(time,state,control,out,grad,args)
 *  \var Diff::scargs
 *  Additional arguments to dynamics
 *  \var Diff::s
 *  RHS of diffusion term of stochastic differential equation
 *  f(time,state,out,args)
 *  \var Diff::sargs
 *  Additional arguments to dynamics
 */
struct Diff
{
    
    int controlled;
    size_t dx;
    size_t du;
    size_t dw;

    int (*sc)(double,const double*,const double*,double*,double*,void*);
    void * scargs;

    int (*s)(double,const double*,double*,void*);
    void * sargs;
};


/*!
  Add controlled drift dynamcis

  \param[in,out] dr    - allocated diff object
  \param[in]     s     - controlled dynamics
  \param[in]     sargs - arguments to dynamics

*/
void diff_add_rhs_controlled(struct Diff * diff,
                             int (*s)(double,const double*,const double*,
                                      double*,double*,void*),
                             void * sargs)
{
    assert (diff != NULL);
    diff->sc = s;
    diff->scargs = sargs;
}

/*!
  Allocate controlled diffusion dynamics

  \param[in] dx     - dimension of state
  \param[in] du     - dimension of control
  \param[in] dw     - dimension of noise
  \param[in] sd     - controlled dynamics
  \param[in] sdargs - arguments to dynamics

  \return diffusion dynamics
*/
struct Diff * 
diff_alloc_controlled(size_t dx, size_t du, size_t dw,
                      int (*sd)(double,const double*,const double*,
                                double*,double*,void*),
                      void * sdargs)
{
    struct Diff * s = malloc(sizeof(struct Diff));
    if (s == NULL){
        fprintf(stderr,"Error allocating memory for Diffusion\n");
        exit(1);
    }

    s->controlled = 1;
    s->dx = dx;
    s->du = du;
    s->dw = dw;

    diff_add_rhs_controlled(s,sd,sdargs);
    s->s = NULL;
    s->sargs = NULL;

    return s;
}

/*!
  Add uncontrolled drift dynamcis

  \param[in,out] dr    - allocated diff object
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
diff_alloc_uncontrolled(size_t dx, size_t dw,
                        int (*sd)(double,const double*,double*,void*),
                        void * sdargs)
{
    struct Diff * s = malloc(sizeof(struct Diff));
    if (s == NULL){
        fprintf(stderr,"Error allocating memory for Diff\n");
        exit(1);
    }
    
    s->controlled = 0;
    s->dx = dx;
    s->dw = dw;
    s->du = 0;

    s->sc = NULL;
    s->scargs = NULL;

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
              const double * u, double * out, double * jac)
{
    int res;
    if (b->controlled == 0){
        /* assert(u == NULL); */
        /* assert(jac == NULL); */
        assert(b->s != NULL);
        res = b->s(time,x,out,b->sargs);
    }
    else{
        assert(b->sc != NULL);
        res = b->sc(time,x,u,out,jac,b->scargs);
    }

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
   Get size of control space
*/
size_t dyn_get_du(const struct Dyn * dyn)
{
    assert (dyn != NULL);
    return drift_get_du(dyn->drift);
}

/*!
   Evaluate dynamics
*/
int dyn_eval(const struct Dyn * dyn,double time,
             const double * x,
             const double * u, double * drift, 
             double * jacdr,
             double * diff, double * jacdiff)
{
    int res = 0;
    if (drift != NULL){
        res = drift_eval(dyn->drift,time,x,u,drift,jacdr);
    }
    if (res != 0){
        return res;
    }
    if (diff != NULL){
        res = diff_eval(dyn->diff,time,x,u,diff,jacdiff);
    }
    return res;
}


