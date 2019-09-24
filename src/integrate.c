// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "integrate.h"

enum IntegratorType {UNSET,FE,RK4,RKF45};

/** \struct Integrator
 *  \brief Integrator
 *  \var Integrator::dx
 *  dimension of state space
 *  \var Integrator::du
 *  dimension of control space
 *  \var Integrator::drift
 *  drift dynamics (time,state,out,jac,args)
 *  \var Integrator::dargs
 *  additional arguments to drift dynamics
 *  \var Integrator::dt
 *  timestep to use
 *  \var Integrator::dtmin
 *  minimum timestep
 *  \var Integrator::dtmax
 *  minimum timestep
 *  \var Integrator::tol
 *  tolerance for adaptive schemes
 *  \var Integrator::space
 *  allocated space
 *  \var Integrator::verbose
 *  verbosity level
 */
struct Integrator
{
    size_t dx;

    enum IntegratorType type;
    int (*drift)(double,const double *,double*,double*,void*);
    void * dargs;

    // controller stuff
    size_t du;
    int (*controlled_dyn)(double, const double *,
                          const double *, double *,
                          double *, void *);
    void * cdargs;
    int (*controller)(double, const double *, double *, void *);
    void * cargs;
    
    // for fixed step
    double dt;
    
    // these are necessary for adaptive schemes
    double dtmin;
    double dtmax;
    double tol;

    double * space;

    int verbose;
};

/*!
  Create an integrator

  \param[in] dx    - dimension of state
  \param[in] b     - drift dynamics(time,x,out,jac,bargs)
  \param[in] bargs - arguments to drift dynamics

  \return integrator
*/
struct Integrator *
integrator_create(size_t dx,
                  int (*b)(double,const double*,double*,double*,void*),
                  void * bargs)
{
    struct Integrator * i = malloc(sizeof(struct Integrator));
    if (i == NULL){
        fprintf(stderr,"Error allocating memory for Integrator\n");
        exit(1);
    }
    i->dx = dx;
    i->type = UNSET;
    i->drift = b;
    i->dargs = bargs;
    i->space = NULL;

    i->du = 0;
    i->controlled_dyn = NULL;
    i->cdargs = NULL;
    i->controller = NULL;
    i->cargs = NULL;
    
    i->dt = 0;
    i->dtmin = 0;
    i->dtmax = 0;
    i->tol = 0.0;

    i->verbose = 0;
    return i;
}

int bcon(double time, const double * x, double * out, double * jac, void * arg)
{
    assert (jac==NULL);
    struct Integrator * i = arg;
    double * u = calloc(i->du,sizeof(double));
    if (u == NULL){
        return 1;
    }
    int res = i->controller(time,x,u,i->cargs);
    if (res != 0 ){
        printf("Controller returned bad value\n");
        free(u); u= NULL;
        return res;
    }
    res = i->controlled_dyn(time,x,u,out,NULL,i->cdargs);
    free(u); u = NULL;
    return res;
}

/*!
  Create an integrator

  \param[in] dx         - dimension of state
  \param[in] du         - dimension of control
  \param[in] b          - drift dynamics(time,x,u,out,jac,bargs)
  \param[in] bargs      - arguments to drift dynamics
  \param[in] controller - controller(time,x,out,cargs) 
  \param[in] cargs      - control arguments

  \return integrator
*/
struct Integrator *
integrator_create_controlled(
    size_t dx, size_t du,
    int (*b)(double, const double *, const double *,double*,double *, void *),
    void * bargs,
    int (*controller)(double, const double *, double *, void *),
    void * cargs)
{

    struct Integrator * i = malloc(sizeof(struct Integrator));
    if (i == NULL){
        fprintf(stderr,"Error allocating memory for Integrator\n");
        exit(1);
    }
    i->dx = dx;
    i->type = UNSET;
    i->space = NULL;

    i->du = du;
    i->controlled_dyn = b;
    i->cdargs = bargs;
    i->controller = controller;
    i->cargs = cargs;

    i->dt = 0;
    i->dtmin = 0;
    i->dtmax = 0;
    i->tol = 0.0;

    i->drift = bcon;
    i->dargs = i; // this is extremly dangerous!!!

    return i;
}

/*!
Evaluate the controller
*/
int integrator_eval_controller(struct Integrator * i,
                               double time, const double * x, double * u)
{
    int res = i->controller(time,x,u,i->cargs);
    return res;
}

/*!
Free memory for integrator
*/
void integrator_destroy(struct Integrator * i)
{
    if (i != NULL){
        free(i->space); i->space = NULL;
        free(i); i= NULL;
    }
}

/*!
  Set type of integrator

  \param[in,out] i    - integrator
  \param[in]     type - forward-euler, rk4, rk45
  
*/
void integrator_set_type(struct Integrator * i, char * type)
{

    if (strcmp(type,"forward-euler") == 0){
        i->type = FE;
        i->space = malloc(i->dx*sizeof(double));
    }
    else if (strcmp(type,"rk4") == 0){
        i->type = RK4;
        i->space = malloc(4*i->dx*sizeof(double));
    }
    else if (strcmp(type,"rkf45") == 0){
        i->type = RKF45;
        i->space = malloc(6*i->dx*sizeof(double));
    }
    else{
        fprintf(stderr,"Integrator type %s is not known\n",type);
    }

    if (i->space == NULL){
        fprintf(stderr, "Memory error allocating space for integrator of type %s\n",type);
        exit(1);
    }
}

/*!
  Set time step of fixed step integrators 

  \param[in,out] i  - integrator
  \param[in]     dt - time_step
*/
void integrator_set_dt(struct Integrator * i, double dt)
{
    assert (i!= NULL);
    i->dt = dt;
}

/*!
  Set minimum and maximum time step for variable step integrators

  \param[in,out] i     - integrator
  \param[in]     dtmin - minimum timestep
  \param[in]     dtmax - maximum timestep
  \param[in]     tol   - tolerance
*/
void integrator_set_adaptive_opts(struct Integrator * i, double dtmin, double dtmax, double tol)
{
    assert (i!= NULL);
    i->dtmin = dtmin;
    i->dtmax = dtmax;
    i->tol = tol;
}

/*!
  Set verbosity level
*/
void integrator_set_verbose(struct Integrator * i, int verbose)
{
    assert (i!= NULL);
    i->verbose = verbose;
}

/*!
  Set the drift arguments

  \param[in,out]   - integrator
  \param[in] bargs - arguments to drift dynamics

*/
void integrator_set_dargs(struct Integrator * ode, void * bargs)
{
    assert (ode != NULL);
    ode->dargs = bargs;
}

/*!
  Fourth order runge kutta scheme
*/
int integrator_step_rk4(struct Integrator * i, double time, double * x)
{
        

    int res = i->drift(time,x,i->space,NULL,i->dargs);
    if (res != 0){
        return res;
    }

    //k2
    for (size_t ii = 0; ii < i->dx; ii++){
        x[ii] += i->dt/2.0 * i->space[ii];
    }

    res = i->drift(time+i->dt/2.0,x,i->space+i->dx,NULL,i->dargs);
    if (res != 0){
        return res;
    }

    //k3
    for (size_t ii = 0; ii < i->dx; ii++){
        x[ii] -= i->dt/2.0 * i->space[ii];
        x[ii] += i->dt/2.0 * i->space[ii+i->dx];
    }

    res = i->drift(time+i->dt/2.0,x,i->space+2*i->dx,NULL,i->dargs);
    if (res != 0){
        return res;
    }

    //k4
    for (size_t ii = 0; ii < i->dx; ii++){
        x[ii] -= i->dt/2.0 * i->space[ii+i->dx];
        x[ii] += i->dt * i->space[ii+2*i->dx];
    }
    res = i->drift(time+i->dt/2.0,x,i->space+3*i->dx,NULL,i->dargs);
    if (res != 0){
        return res;
    }

    // final
    for (size_t ii = 0; ii < i->dx; ii++){
        x[ii] -= i->dt * i->space[ii+2*i->dx];
        x[ii] += i->dt/6.0 * ( i->space[ii] + 2.0*i->space[i->dx+ii] +
                                  2.0*i->space[2*i->dx+ii] + i->space[3*i->dx+ii]);
    }

    return res;
}


/*!
  Fourth order / fifth-order runge kutta scheme (RKF45)
*/
double integrator_step_rkf45(struct Integrator * i, double time, double * x)
{
    double c[6] = {0.0, 0.25, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0};
    double a1 = 1.0/4.0;
    double a2[2] = {3.0/32.0, 9.0/32.0};
    double a3[3] = {1932.0/2197.0,-7200.0/2197.0,7296.0/2197.0};
    double a4[4] = {439.0/216.0,-8.0,3680.0/513.0,-845.0/4104.0};
    double a5[5] = {-8.0/27.0,2.0,-3544.0/2565.0,1859.0/4104.0,-11.0/40.0};
//    double fourth[6] = {25.0/216.0,0.0,1408.0/2565.0,2197.0/4104.0,-1.0/5.0,0.0};
    double fifth[6] = {16.0/135.0,0.0,6656.0/12825.0,28561.0/56430.0,-9.0/50.0,2.0/55.0};
    double * xs = malloc(i->dx * sizeof(double));
    int converged = 0;
    double next_dt;
    do
    {
        //k1
        int res = i->drift(time+c[0]*i->dt,x,i->space,NULL,i->dargs);
        if (res != 0){
            return res;
        }
        
        //k2
        for (size_t ii = 0; ii < i->dx; ii++){
            i->space[ii] *= i->dt;
            xs[ii] = x[ii] + a1 * i->space[ii];

        }
        res = i->drift(time+c[1]*i->dt,xs,i->space+i->dx,NULL,i->dargs);
        if (res != 0){
            return res;
        }

        //k3
        for (size_t ii = 0; ii < i->dx; ii++){
            xs[ii] =  x[ii] + a2[0]* i->space[ii];
            i->space[ii + i->dx] *= i->dt;
            xs[ii] += a2[1] * i->space[ii+i->dx];

        }
        res = i->drift(time+c[2]*i->dt,xs,i->space+2*i->dx,NULL,i->dargs);
        if (res != 0){
            return res;
        }

        //k4
        for (size_t ii = 0; ii < i->dx; ii++){
            xs[ii] =  x[ii] + a3[0] * i->space[ii];
            xs[ii] += a3[1] * i->space[ii+i->dx];
            i->space[ii + 2*i->dx] *= i->dt;
            xs[ii] += a3[2] * i->space[ii+2*i->dx];
        }
        res = i->drift(time+c[3]*i->dt,xs,i->space+3*i->dx,NULL,i->dargs);
        if (res != 0){
            return res;
        }

        //k5
        for (size_t ii = 0; ii < i->dx; ii++){
            xs[ii] =  x[ii] + a4[0] * i->space[ii];
            xs[ii] += a4[1] * i->space[ii+i->dx];
            xs[ii] += a4[2] * i->space[ii+2*i->dx];
            i->space[ii + 3*i->dx] *= i->dt;
            xs[ii] += a4[3] * i->space[ii+3*i->dx];
        }
        res = i->drift(time+c[4]*i->dt,xs,i->space+4*i->dx,NULL,i->dargs);
        if (res != 0){
            return res;
        }

        //k6
        for (size_t ii = 0; ii < i->dx; ii++){
            xs[ii] =  x[ii] + a5[0] * i->space[ii];
            xs[ii] += a5[1] * i->space[ii+i->dx];
            xs[ii] += a5[2] * i->space[ii+2*i->dx];
            xs[ii] += a5[3] * i->space[ii+3*i->dx];
            i->space[ii + 4*i->dx] *= i->dt;
            xs[ii] += a5[4] * i->space[ii+4*i->dx];
        }
        res = i->drift(time+c[5]*i->dt,xs,i->space+5*i->dx,NULL,i->dargs);
        if (res != 0){
            return res;
        }

        double abserr = 0.0;
        for (size_t ii = 0; ii < i->dx; ii++){
            i->space[ii + 5*i->dx] *= i->dt;
            
            /* printf("ii = %zu\n",ii); */
            /* printf("k1=%G\n",i->space[ii]); */
            /* printf("k2=%G\n",i->space[ii+i->dx]); */
            /* printf("k3=%G\n",i->space[ii+2*i->dx]); */
            /* printf("k4=%G\n",i->space[ii+3*i->dx]); */
            /* printf("k5=%G\n",i->space[ii+4*i->dx]); */
            /* printf("k6=%G\n",i->space[ii+5*i->dx]); */
            double err = 1.0/360 * i->space[ii];
//            printf("err = %G\n",err);
            err -= 128.0/4275.0 * i->space[ii + 2*i->dx];
//            printf("err = %G\n",err);
            err -= 2197.0/75240.0 * i->space[ii + 3*i->dx];
//            printf("err = %G\n",err);
            err += 1.0/50.0 * i->space[ii+4*i->dx];
//            printf("err = %G\n",err);
            err += 2.0/55.0 * i->space[ii+5*i->dx];
//            printf("err = %G\n",err);
            err /= i->dt;
//            printf("err = %G\n",err);
            abserr += fabs(err);
//            for (size_t jj = 0; jj < 6; jj++){
                //y[ii] += fourth[jj]*i->space[ii + jj*i->dx];
//                z[ii] += fifth[jj]*i->space[ii + jj*i->dx];
                //          }
//            err[ii] = z[ii]- y[ii];
//            abserr += pow(err[ii],2);
        }

        if (isnan(abserr)){
            printf("abserr is NAN\n");
            exit(1);
        }
        double q = pow(i->tol/(2.0 * abserr),0.25);

        next_dt = q * i->dt;
        if (next_dt > i->dtmax){
            next_dt = i->dtmax;
        }
        else if (next_dt < i->dtmin){
            next_dt = i->dtmin;
        }

        if (i->verbose > 1){

            fprintf(stdout,"abserr=%G, tol=%G\n",abserr,i->tol);
            fprintf(stdout,"Current dt=%G, next_dt = %G\n",i->dt,next_dt);
            fprintf(stdout,"q = %G\n",q);
        }
        
        if (abserr < i->tol){
            converged = 1;
        }
        else if (fabs(next_dt - i->dt) < 1e-14){
            converged = 1;
        }
        else{
            i->dt = next_dt;
        }

    }while (converged == 0);

    for (size_t ii = 0; ii < i->dx; ii++){
        for (size_t jj = 0; jj < 6; jj++){
            x[ii] += fifth[jj]*i->space[ii + jj*i->dx];
        }
    }
    
    free(xs); xs = NULL;
    return next_dt;
}

/*!
  Integrate from *start_time* to *end_time*

  \param[in] i          - integrator
  \param[in] start_time - initial time
  \param[in] end_time   - ending time (will integrate until reach past this time)
  \param[in] start      - starting state
  \param[in] end        - ending state
*/
void integrator_step(struct Integrator * i,
                     double start_time,
                     double end_time,
                     const double * start, double * end)
{

    assert (i != NULL);
    assert (i->type != UNSET);
    
    double time = start_time;
    if (i->type == FE){
        if (i->dt <= 0){
            fprintf(stderr, "Must set timestep before integrating\n");
            exit(1);
        }
        for (size_t ii = 0; ii < i->dx; ii++){
            end[ii] = start[ii];
        }
        while (time < end_time){
            int res = i->drift(time,end,i->space,NULL,i->dargs);
            assert (res == 0);
            for (size_t ii = 0; ii < i->dx; ii++){
                end[ii] = end[ii] + i->dt * i->space[ii];
            }
            time = time + i->dt;
        }
    }
    else if (i->type == RK4){
         if (i->dt <= 0){
            fprintf(stderr, "Must set timestep before integrating\n");
            exit(1);
        }
        for (size_t ii = 0; ii < i->dx; ii++){
            end[ii] = start[ii];
        }
        while (time < end_time){
            int res = integrator_step_rk4(i,time,end);
            assert (res == 0);
            time = time + i->dt;
        }
    }
    else if (i->type == RKF45){
        if (i->dtmin <= 0){
            fprintf(stderr, "Must set smallest and largest timesteps and tolerance for adaptive scheme\n");
            exit(1);
        }
        i->dt = i->dtmax;
        for (size_t ii = 0; ii < i->dx; ii++){
            end[ii] = start[ii];
        }
        while (time < end_time){
            
            if (i->verbose > 0){
                fprintf(stdout,"Time=%G\n",time);
            }
            double dtmax_store = i->dtmax;
            if ((end_time-time) < i->dtmax){
//                printf("adjust final time\n");
                i->dtmax = end_time-time;
                if (i->dt > i->dtmax){
                    i->dt = i->dtmax;
                }
            }
            double next_dt = integrator_step_rkf45(i,time,end);
            time = time + i->dt;
            i->dt = next_dt;
            i->dtmax = dtmax_store;
        }
        if (i->verbose > 0){
            fprintf(stdout,"Finished Time=%G\n",time);
        }
    }
}

/*!
   Take a step of euler_maruyama
*/
/* int  */
/* euler_maruyama_step(double time, const double * x, const double * noise, */
/*                     double * nextx, double dt, struct Dyn * dyn,  */
/*                     double * drift, double * diff) */
/* { */

/*     size_t d = dyn_get_dx(dyn); */
/*     size_t dw = dyn_get_dw(dyn); */
/*     int res = dyn_eval(dyn,time,x,drift,NULL,diff); */

/*     double sqrtdt = sqrt(dt); */
/*     for (size_t ii = 0; ii < d; ii++ ){ */
/*         nextx[ii] = x[ii] + dt * drift[ii]; */
/*         for (size_t jj = 0; jj < dw; jj++){ */
/*             nextx[ii] += sqrtdt*diff[jj*d+ii]*noise[jj]; */
/*             if (isinf(nextx[ii])){ */
/*                 fprintf(stderr,"STOP!\n"); */
/*                 printf("noise \n"); */
/*                 for (size_t kk = 0; kk < dw; kk++){ */
/*                     printf("%G ",noise[kk]); */
/*                 } */
/*                 printf("\n"); */
/*                 printf("now diff\n"); */
/*                 for (size_t kk = 0; kk < dw; kk++){ */
/*                     printf("%G ",diff[kk*d+ii]); */
/*                 } */
/*                 printf("\n"); */
/*                 printf("drift is %G\n",drift[ii]); */
/*                 printf("now x\n"); */
/*                 for (size_t kk = 0; kk < d; kk++){ */
/*                     printf("%G ",x[kk]); */
/*                 } */
/*                 printf("\n"); */
/*                 exit(1); */
/*             } */
/*         } */
            
/*     } */
/*     return res; */
/* } */



