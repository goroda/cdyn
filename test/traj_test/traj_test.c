// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice, 
//    this list of conditions and the following disclaimer in the documentation 
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its contributors 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//Code


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "CuTest.h"
#include "integrate.h"
#include "simulate.h"

static size_t dx = 3;
static size_t du = 0;

////////////////////////////////////////////////////////////////////

void Test_trajectory_alloc(CuTest * tc)
{

    fprintf(stdout,"Testing allocation and additions to trajectory\n");

    double t1 = 0.5;
    double pt1[3] = {0.5, 0.2, 0.1};

    double t2 = 1.5;
    double pt2[3] = {-0.5, 0.5, 2.0};

    double t3 = 2.5;
    double pt3[3] = {-0.5, 0.2, 3.2};

    int res;
    struct Trajectory * traj = NULL;
    res = trajectory_add(&traj,dx,du,t1,pt1,NULL);
    CuAssertIntEquals(tc,0,res);
    res = trajectory_add(&traj,dx,du,t2,pt2,NULL);
    CuAssertIntEquals(tc,0,res);
    res = trajectory_add(&traj,dx,du,t3,pt3,NULL);
    CuAssertIntEquals(tc,0,res);

    double last_time = trajectory_get_last_time(traj);
    CuAssertDblEquals(tc,t3,last_time,1e-20);

    double * lp = trajectory_get_last_state(traj);
    double * lu = trajectory_get_last_control(traj);
    for (size_t ii = 0; ii < dx; ii++){
        CuAssertDblEquals(tc,pt3[ii],lp[ii],1e-20);
    }
    CuAssertIntEquals(tc,1,lu==NULL);

//    trajectory_print(traj,stdout,2);
    trajectory_free(traj); traj = NULL;
}

int b(double time, const double * x, double * out, double * jac, void * arg)
{
    assert (jac == NULL);
    (void)(arg);
    out[0] = time*x[0];
    out[1] = x[1];
    out[2] = x[2]*x[2];
    return 0;
}

void Test_trajectory_step(CuTest * tc)
{

    fprintf(stdout,"Testing trajectory_step\n");

    struct Integrator * ode = integrator_create(dx,b,NULL);
    double dtmin = 1e-3;
    double dtmax = 0.25;
    double tol = 1e-8;
    integrator_set_type(ode,"rkf45");
    integrator_set_adaptive_opts(ode,dtmin,dtmax,tol);
    integrator_set_verbose(ode,0);

    double t1 = 0.0;
    double pt1[3] = {0.5, 0.2, 0.1};

    int res;
    struct Trajectory * traj = NULL;
    res = trajectory_add(&traj,dx,0,t1,pt1,NULL);
    CuAssertIntEquals(tc,0,res);

    double dt = 1e-1;
    size_t nsteps = 10;
    for (size_t ii = 0; ii < nsteps; ii++){
        res = trajectory_step(traj,ode,dt);
        CuAssertIntEquals(tc,0,res);
    }

    double final_time = dt*nsteps;;
    double sol1 = pt1[0] *exp(pow(final_time,2)/2.0);
    double sol2 = pt1[1] *exp(final_time);
    double sol3 = pt1[2]/(1.0-pt1[2]*final_time);
    
    double * sol = trajectory_get_last_state(traj);
    
    CuAssertDblEquals(tc,sol1,sol[0],1e-2);
    CuAssertDblEquals(tc,sol2,sol[1],1e-2);
    CuAssertDblEquals(tc,sol3,sol[2],1e-2);
    
//    trajectory_print(traj,stdout,2);
    trajectory_free(traj); traj = NULL;
    integrator_destroy(ode); ode = NULL;
}

CuSuite * TrajGetSuite(){

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_trajectory_alloc);
    SUITE_ADD_TEST(suite, Test_trajectory_step);
    return suite;
}

void RunAllTests(void) {
    
    printf("Running Test Suite: control\n");

    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * traj = TrajGetSuite();

    CuSuiteAddSuite(suite,traj);
        
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);

    CuSuiteDelete(traj);
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}
