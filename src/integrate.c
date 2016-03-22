// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "integrate.h"
#include "dynamics.h"
//#include "simulate.h"

/*!
   Take a step of forward euler
*/
int euler_step(double time, const double * x, double * nextx,
               double dt,struct Dyn * dyn, double * driftv)
{

    int res = dyn_eval(dyn,time,x,NULL,driftv,NULL,NULL,NULL);
       
    size_t dx = dyn_get_dx(dyn);
    for (size_t ii = 0; ii < dx; ii++){
        nextx[ii] = x[ii] + dt * driftv[ii];
    }
    
    return res;
}

int 
euler_maruyama_step(double time, const double * x, const double * noise,
                    double * nextx, double dt, struct Dyn * dyn, 
                    double * drift, double * diff)
{

    size_t d = dyn_get_dx(dyn);
    size_t dw = dyn_get_dw(dyn);
    int res = dyn_eval(dyn,time,x,NULL,drift,NULL,diff,NULL);

    double sqrtdt = sqrt(dt);
    for (size_t ii = 0; ii < d; ii++ ){
        nextx[ii] = x[ii] + dt * drift[ii];
        for (size_t jj = 0; jj < dw; jj++){
            nextx[ii] += sqrtdt*diff[jj*d+ii]*noise[jj];
            if (isinf(nextx[ii])){
                fprintf(stderr,"STOP!\n");
                printf("noise \n");
                for (size_t kk = 0; kk < dw; kk++){
                    printf("%G ",noise[kk]);
                }
                printf("\n");
                printf("now diff\n");
                for (size_t kk = 0; kk < dw; kk++){
                    printf("%G ",diff[kk*d+ii]);
                }
                printf("\n");
                printf("drift is %G\n",drift[ii]);
                printf("now x\n");
                for (size_t kk = 0; kk < d; kk++){
                    printf("%G ",x[kk]);
                }
                printf("\n");
                exit(1);
            }
        }
            
    }
    return res;
}



