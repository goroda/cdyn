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
    \param[in]     ode    - integrator
    \param[in]     dt     - time step
**************************************************************/
int trajectory_step(struct Trajectory * traj,             
                    struct Integrator * ode,
                    double dt)
                    /* struct Dyn * dyn, */
                    /* double dt, char * method, */
                    /* double * space, */
                    /* void * noise, void * args) */
{
    if (traj == NULL){
        fprintf(stderr,"Warning: cannot advance trajectory starting\n");
        fprintf(stderr,"         from NULL Trajectory\n");
        return 1;
    }

    assert (dt > 0);

    double curr_time = trajectory_get_last_time(traj);
    double * x = trajectory_get_last_state(traj);
    size_t dx = traj->dx;
    double * newx = malloc(dx*sizeof(double));
    assert (newx != NULL);

    integrator_step(ode,curr_time,curr_time+dt,x,newx);
    double new_time = curr_time + dt;
    int res = trajectory_add(&traj,dx,0,new_time,newx,NULL);
    free(newx); newx = NULL;;

    return res;

}
