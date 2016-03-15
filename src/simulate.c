// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "dynamics.h"
#include "integrate.h"

/** \struct Trajectory
 *  \brief Linked List defining a trajectory
 *  \var Trajectory::dx
 *  dimension of state space
 *  \var Trajectory::du
 *  dimension of control
 *  \var Trajectory::time
 *  time
 *  \var Trajectory::state
 *  current state
 *  \var Trajectory::control
 *  current control
 *  \var Trajectory::next
 *  Pointer to next element
 */
struct Trajectory
{
    size_t dx;
    size_t du;
    double time;
    double * state;
    double * control;
    struct Trajectory * next;
};

/**********************************************************//**
    Allocate Trajectory
**************************************************************/
struct Trajectory * trajectory_init(size_t dx, size_t du, double time, 
                                    const double * state, 
                                    const double * control)
{

    struct Trajectory * t = malloc(sizeof(struct Trajectory ));
    if (t == NULL){
        fprintf(stderr, "Cannot allocater ajectory\n");
        exit(1);
    }
    
    t->dx = dx;
    t->du = du;
    t->time = time;
    t->state = malloc(dx*sizeof(double));
    if (t->state == NULL){
        fprintf(stderr, "Cannot allocate state in Trajectory\n");
        exit(1);
    }
    memmove(t->state,state,dx*sizeof(double));
    if (du != 0){
        t->du = du;
        t->control = malloc(du*sizeof(double));
        if (t->control == NULL){
            fprintf(stderr, "Cannot allocate control in Trajectory\n");
            exit(1);
        }
        memmove(t->control,control,du*sizeof(double));
    }
    else{
        t->du = 0;
        t->control = NULL;
    }
    t->next = NULL;
    return t;
}

/**********************************************************//**
    Free Trajectory
**************************************************************/
void trajectory_free(struct Trajectory * traj)
{
    if (traj != NULL){
        trajectory_free(traj->next); traj->next = NULL;
        free(traj->state); traj->state = NULL;
        free(traj->control); traj->control = NULL;
        free(traj); traj = NULL;
    }
}

/**********************************************************//**
    Add a state and control to the trajectory

    \param[in,out] traj - trajectory
    \param[in]     dx   - dimension of state
    \param[in]     du   - dimension of control
    \param[in]     t    - time
    \param[in]     x    - state
    \param[in]     u    - control
**************************************************************/
int trajectory_add(struct Trajectory ** traj, size_t dx, size_t du,
                   double t,const double * x,const double * u)
{
    if (*traj == NULL){
        (*traj) = trajectory_init(dx,du,t,x,u);
    }
    else{
        trajectory_add(&((*traj)->next),dx,du,t,x,u);
    }
    return 0;
}

/**********************************************************//**
    Print the trajectory

    \param[in,out] traj - trajectory
    \param[in]     fp   - stream to print to
    \param[in]     prec - precision with which to print
**************************************************************/
void trajectory_print(const struct Trajectory * traj,
                      FILE *fp,
                      int prec)
{
    
    if (traj == NULL){
//        fprintf(fp,"NULL\n");
        fprintf(fp,"\n");
    }
    else{
        fprintf(fp,"%3.*G ",prec,traj->time);
        for (size_t ii = 0; ii < traj->dx; ii++){
            fprintf(fp,"%3.*G ",prec,traj->state[ii]);
        }

        for (size_t ii = 0; ii < traj->du; ii++){
            fprintf(fp,"%3.*G ",prec,traj->control[ii]);
        }
        fprintf(fp,"\n");
        trajectory_print(traj->next,fp,prec);
    }
}

/**********************************************************//**
    Get the last time                                                           
**************************************************************/
double trajectory_get_last_time(const struct Trajectory * traj)
{
    if (traj == NULL){
        fprintf(stderr,"Warning: Cannot get last time, Trajectory is NULL\n");
        return 0.0;
    }
    while (traj->next != NULL){
        traj = traj->next;
    }
    return traj->time;
}

/**********************************************************//**
    Get the last state
**************************************************************/
double * trajectory_get_last_state(const struct Trajectory * traj)
{
    if (traj == NULL){
        fprintf(stderr,"Warning: Cannot get last time, Trajectory is NULL\n");
        return NULL;
    }
    while (traj->next != NULL){
        traj = traj->next;
    }
    return traj->state;
}

/**********************************************************//**
    Get the last control
**************************************************************/
double * trajectory_get_last_control(const struct Trajectory * traj)
{
    if (traj == NULL){
        fprintf(stderr,"Warning: Cannot get last time, Trajectory is NULL\n");
        return NULL;
    }
    while (traj->next != NULL){
        traj = traj->next;
    }
    return traj->control;
}

/**********************************************************//**
    Take a step of a trajectory

    \param[in,out] traj   - trajectory
    \param[in]     dyn    - dynamics
    \param[in]     dt     - time step
                            Used when method=
                            "euler" or 
                            "euler-maruyama"
    \param[in]     method - integration algorithm
    \param[in]     space  - free space to use for computation
    \param[in]     noise  - noise arguments
    \param[in]     args   - additional arguments

    \note
    Let d denote dimension of state 
    Let dw denote dimension of noise 
    Method "euler" - space (d), 
                     noise NULL, 
                     args NULL  
    Method "euler-maruyama" - space(d + d*dw), 
                              noise (dw) gaussian samples,  
                              args NULL 
**************************************************************/
int trajectory_step(struct Trajectory * traj,
                    struct Dyn * dyn,
                    double dt, char * method,
                    double * space,
                    void * noise, void * args)
{
    (void)(args);
    if (traj == NULL){
        fprintf(stderr,"Warning: cannot advance trajectory starting\n");
        fprintf(stderr,"         from NULL Trajectory\n");
        return 1;
    }
    if (dyn == NULL){
        fprintf(stderr,"Warning: cannot advance trajectory starting\n");
        fprintf(stderr,"         from NULL Dyn\n");
        return 1;
    }
    assert (dt > 0);

    
    double curr_time = trajectory_get_last_time(traj);
    double * x = trajectory_get_last_state(traj);
//    double * u = trajectory_get_last_control(traj);

    size_t dx = dyn_get_dx(dyn);
    double * newx = malloc(dx*sizeof(double));
    assert (newx != NULL);

    int res;
    if (strcmp(method,"euler") == 0){
        res = euler_step(curr_time,x,newx,dt,dyn,space);
    }
    else if (strcmp(method,"euler-maruyama") == 0){
        res = euler_maruyama_step(curr_time,x,noise,newx,dt,
                                dyn,space,space+dx);
    }
    else{
        return 1;
    }

    if (res != 0){
        free(newx); newx = NULL;
        return res;
    }

    size_t du = dyn_get_du(dyn);
    double new_time = curr_time + dt;
    res = trajectory_add(&traj,dx,du,new_time,newx,NULL);
    free(newx); newx = NULL;;

    return res;

}
