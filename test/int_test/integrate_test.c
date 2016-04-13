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
//#include "simulate.h"

static size_t dx = 3;

////////////////////////////////////////////////////////////////////
int b(double time, const double * x, double * out, double * jac, void * arg)
{
    assert (jac == NULL);
    (void)(arg);
    out[0] = time*x[0];
    out[1] = x[1];
    out[2] = x[2]*x[2];
    return 0;
}

void Test_integrator_forward_euler(CuTest * tc)
{

    fprintf(stdout,"Testing forward_euler\n");

    struct Integrator * ode = integrator_create(dx,b,NULL);
    integrator_set_type(ode,"forward-euler");
    integrator_set_dt(ode,1e-3);
    
    double start_time = 0.0;
    double end_time = 2.0;
    double pt1[3] = {0.5, 0.2, 0.1};
    double sol[3];

    integrator_step(ode,start_time,end_time,pt1,sol);

    double final_time = end_time;;
    double sol1 = pt1[0] *exp(pow(final_time,2)/2.0);
    double sol2 = pt1[1] *exp(final_time);
    double sol3 = pt1[2]/(1.0-pt1[2]*final_time);
    
    CuAssertDblEquals(tc,sol1,sol[0],1e-2);
    CuAssertDblEquals(tc,sol2,sol[1],1e-2);
    CuAssertDblEquals(tc,sol3,sol[2],1e-2);

    integrator_destroy(ode);
}

void Test_integrator_rk4(CuTest * tc)
{

    fprintf(stdout,"Testing rk4\n");

    struct Integrator * ode = integrator_create(dx,b,NULL);
    integrator_set_type(ode,"rk4");
    integrator_set_dt(ode,1e-3);
    
    double start_time = 0.0;
    double end_time = 2.0;
    double pt1[3] = {0.5, 0.2, 0.1};
    double sol[3];

    integrator_step(ode,start_time,end_time,pt1,sol);

    double final_time = end_time;;
    double sol1 = pt1[0] *exp(pow(final_time,2)/2.0);
    double sol2 = pt1[1] *exp(final_time);
    double sol3 = pt1[2]/(1.0-pt1[2]*final_time);
    
    CuAssertDblEquals(tc,sol1,sol[0],1e-2);
    CuAssertDblEquals(tc,sol2,sol[1],1e-2);
    CuAssertDblEquals(tc,sol3,sol[2],1e-2);

    integrator_destroy(ode);
}

int b45(double time, const double * x, double * out, double * jac, void * arg)
{
    assert (jac == NULL);
    (void)(arg);
    out[0] = -3*time*pow(x[0],2) + 1.0/(1.0 + pow(time,3));
//    out[1] = x[1];
//    out[2] = x[2]*x[2];
    return 0;
}

void Test_integrator_rkf45(CuTest * tc)
{

    fprintf(stdout,"Testing rkf45\n");

//    double dtmin = 1e-2;
    double dtmin = 1e-3;
    double dtmax = 0.25;
    double tol = 5e-8;

//    struct Integrator * ode = integrator_create(1,b45,NULL);
    struct Integrator * ode = integrator_create(dx,b,NULL);
    integrator_set_type(ode,"rkf45");
    integrator_set_adaptive_opts(ode,dtmin,dtmax,tol);
    integrator_set_verbose(ode,0);
    
    /* double start_time = 0.0; */
    /* double end_time = 3e-1; */
    /* double pt1[1] = {0.0}; */
    /* double sol[1]; */
    
    double start_time = 0.0;
    double end_time = 2.0;
    double pt1[3] = {0.5, 0.2, 0.1};
    double sol[3];

    integrator_step(ode,start_time,end_time,pt1,sol);

    double final_time = end_time;;
    double sol1 = pt1[0] *exp(pow(final_time,2)/2.0);
    double sol2 = pt1[1] *exp(final_time);
    double sol3 = pt1[2]/(1.0-pt1[2]*final_time);
    
    CuAssertDblEquals(tc,sol1,sol[0],1e-7);
    CuAssertDblEquals(tc,sol2,sol[1],1e-7);
    CuAssertDblEquals(tc,sol3,sol[2],1e-7);

    integrator_destroy(ode);
}

int b45_cont(double time, const double * x,
             const double * u, double * out, double * jac, void * arg)
{
    assert (jac == NULL);
    (void)(arg);
    out[0] = time*x[0];
    out[1] = x[1];
    out[2] = u[0];
    return 0;
}

int controller(double time, const double * x, double * u, void * arg)
{
    (void)(time);
    (void)(arg);
    u[0] = pow(x[2],2);
    return 0;
}

void Test_integrator_rkf45_controlled(CuTest * tc)
{

    fprintf(stdout,"Testing rkf45 with control\n");

//    double dtmin = 1e-2;
    double dtmin = 1e-3;
    double dtmax = 0.25;
    double tol = 5e-8;

//    struct Integrator * ode = integrator_create(1,b45,NULL);
    struct Integrator * ode =
        integrator_create_controlled(dx,1,b45_cont,NULL,controller,NULL);
    integrator_set_type(ode,"rkf45");
    integrator_set_adaptive_opts(ode,dtmin,dtmax,tol);
    integrator_set_verbose(ode,0);
    
    /* double start_time = 0.0; */
    /* double end_time = 3e-1; */
    /* double pt1[1] = {0.0}; */
    /* double sol[1]; */
    
    double start_time = 0.0;
    double end_time = 2.0;
    double pt1[3] = {0.5, 0.2, 0.1};
    double sol[3];

    integrator_step(ode,start_time,end_time,pt1,sol);

    double final_time = end_time;;
    double sol1 = pt1[0] *exp(pow(final_time,2)/2.0);
    double sol2 = pt1[1] *exp(final_time);
    double sol3 = pt1[2]/(1.0-pt1[2]*final_time);
    
    CuAssertDblEquals(tc,sol1,sol[0],1e-7);
    CuAssertDblEquals(tc,sol2,sol[1],1e-7);
    CuAssertDblEquals(tc,sol3,sol[2],1e-7);

    integrator_destroy(ode);
}

CuSuite * IntegratorGetSuite(){

    CuSuite * suite = CuSuiteNew();
//    SUITE_ADD_TEST(suite, Test_trajectory_alloc);
    SUITE_ADD_TEST(suite, Test_integrator_forward_euler);
    SUITE_ADD_TEST(suite, Test_integrator_rk4);
    SUITE_ADD_TEST(suite, Test_integrator_rkf45);
    SUITE_ADD_TEST(suite, Test_integrator_rkf45_controlled);
    return suite;
}

void RunAllTests(void) {
    
    printf("Running Test Suite: control\n");

    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * ode = IntegratorGetSuite();

    CuSuiteAddSuite(suite,ode);
        
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);

    CuSuiteDelete(ode);
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}
