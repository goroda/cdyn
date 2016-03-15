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

#include "CuTest.h"
#include "dynamics.h"

static size_t dx = 3;
static size_t du = 2;
static size_t dw = 2;
int bc(double time, const double *x, const double *u,
       double *out,double * grad, void * arg)
{
    double coeff = *(double *)arg;
    out[0] = exp(-time *10)*x[1]*x[2] + u[1];
    out[1] = coeff + sin(x[1] + x[2]);
    out[2] = x[1] * x[0] + u[0];

    if (grad != NULL)
    {
        grad[0] = 0; grad[3] = 1;
        grad[1] = 0; grad[4] = 0;
        grad[2] = 1; grad[5] = 0;
    }
    return 0;
}

int b(double time, const double * x, double * out, void * arg)
{
    double coeff = *(double *)arg;
    out[0] = exp(-time *10)*x[1]*x[2];
    out[1] = coeff + sin(x[1] + x[2]);
    out[2] = x[1] * x[0];

    return 0;
}

int sc(double time, const double * x, const double * u,
       double * out, double * grad, void * arg)
{
    double coeff = *(double *)arg;
    
    out[0] = x[0];      out[3] = x[0]*sin(x[1]*x[2]);
    out[1] = time*x[2]; out[4] = coeff+u[1];
    out[2] = u[0];      out[5] = u[0]*u[1];

    if (grad != NULL)
    {
        //du
        grad[0] = 0.0; grad[3] = 0.0;
        grad[1] = 0.0; grad[4] = 0.0;
        grad[2] = 1.0; grad[5] = u[1];

        //du1
        grad[6] = 0.0; grad[9]  = 0.0;
        grad[7] = 0.0; grad[10] = 1.0;
        grad[8] = 0.0; grad[11] = u[0];
    }
    return 0;
}

int s(double time, const double * x,double * out, void * arg)
{
    double coeff = *(double *)arg;
    
    out[0] = x[0];           out[3] = x[0]*sin(x[1]*x[2]);
    out[1] = time*x[2];      out[4] = coeff*x[2];
    out[2] = exp(x[0]+x[1]); out[5] = cos(x[0]*x[1]);

    return 0;
}

////////////////////////////////////////////////////////////////////

void Test_drift_controlled(CuTest * tc)
{
    fprintf(stdout,"Testing controlled drift dynamics\n");
    double coeff = 2.0;
    struct Drift * drift = drift_alloc_controlled(dx,du,bc,&coeff);

    double test_pt[3] = {0.3, 2.0, -0.4};
    double test_u[2] = {5.0, -5.0};
    double time = 2.3;

    double out1[3];
    double out2[3];
    double grad1[6];
    double grad2[6];
    
    int res = drift_eval(drift,time,test_pt,test_u,out1,grad1);
    CuAssertIntEquals(tc,0,res);

    bc(time,test_pt,test_u,out2,grad2,&coeff);
    for (size_t ii = 0; ii < dx; ii++)
    {
        CuAssertDblEquals(tc,out2[ii],out1[ii],1e-20);        
        for (size_t jj = 0; jj < du;jj++){
            size_t ind = jj*dx+ii;
            CuAssertDblEquals(tc,grad2[ind],grad1[ind],1e-20);        
        }
    }
    drift_free(drift);
}

void Test_drift_uncontrolled(CuTest * tc)
{
    fprintf(stdout,"Testing uncontrolled drift dynamics\n");
    double coeff = 2.0;
    struct Drift * drift = drift_alloc_uncontrolled(dx,b,&coeff);

    double test_pt[3] = {0.3, 2.0, -0.4};
    double time = 2.3;

    double out1[3];
    double out2[3];
    
    int res = drift_eval(drift,time,test_pt,NULL,out1,NULL);
    CuAssertIntEquals(tc,0,res);

    b(time,test_pt,out2,&coeff);
    for (size_t ii = 0; ii < dx; ii++)
    {
        CuAssertDblEquals(tc,out2[ii],out1[ii],1e-20);        
    }
    drift_free(drift);
}

void Test_diff_controlled(CuTest * tc)
{
    fprintf(stdout,"Testing controlled diffusion dynamics\n");
    double coeff = 2.0;
    struct Diff * diff = diff_alloc_controlled(dx,du,dw,sc,&coeff);

    double test_pt[3] = {0.3, 2.0, -0.4};
    double test_u[2] = {5.0, -5.0};
    double time = 2.3;

    double out1[6];
    double out2[6];
    double grad1[12];
    double grad2[12];
    
    int res = diff_eval(diff,time,test_pt,test_u,out1,grad1);
    CuAssertIntEquals(tc,0,res);

    sc(time,test_pt,test_u,out2,grad2,&coeff);
    for (size_t ii = 0; ii < dx; ii++)
    {
        for (size_t jj = 0; jj < dw; jj++){
            size_t ind1 = jj*dx+ii;
            CuAssertDblEquals(tc,out2[ind1],out1[ind1],1e-20);        

            for (size_t kk = 0; kk < du;kk++){
                size_t ind2 = kk*dx*dw + ind1;
                CuAssertDblEquals(tc,grad2[ind2],grad1[ind2],1e-20);        
            }
        }
    }
    diff_free(diff);
}

void Test_diff_uncontrolled(CuTest * tc)
{
    fprintf(stdout,"Testing uncontrolled diffusion dynamics\n");
    double coeff = 2.0;
    struct Diff * diff = diff_alloc_uncontrolled(dx,dw,s,&coeff);

    double test_pt[3] = {0.3, 2.0, -0.4};
    double time = 2.3;

    double out1[6];
    double out2[6];
    
    int res = diff_eval(diff,time,test_pt,NULL,out1,NULL);
    CuAssertIntEquals(tc,0,res);

    s(time,test_pt,out2,&coeff);
    for (size_t ii = 0; ii < dx; ii++)
    {
        for (size_t jj = 0; jj < dw; jj++){
            size_t ind1 = jj*dx+ii;
            CuAssertDblEquals(tc,out2[ind1],out1[ind1],1e-20);        
        }
    }
    diff_free(diff);
}

void Test_dyn_alloc(CuTest * tc)
{

    fprintf(stdout,"Testing allocation and evaluation of Dyn\n");
    double coeff = 2.0;
    struct Diff * diff = diff_alloc_uncontrolled(dx,dw,s,&coeff);
    struct Drift * drift = drift_alloc_uncontrolled(dx,b,&coeff);
    struct Dyn * dyn = dyn_alloc(drift,diff);

    double test_pt[3] = {0.3, 2.0, -0.4};
    double time = 2.3;
    double out_drift[3];
    double out_diff[6];

    int res = dyn_eval(dyn,time,test_pt,NULL,out_drift,NULL,out_diff,NULL);
    CuAssertIntEquals(tc,0,res);

    dyn_free_deep(dyn);
}

CuSuite * DynGetSuite(){

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_drift_controlled);
    SUITE_ADD_TEST(suite, Test_drift_uncontrolled);
    SUITE_ADD_TEST(suite, Test_diff_controlled);
    SUITE_ADD_TEST(suite, Test_diff_uncontrolled);
    SUITE_ADD_TEST(suite, Test_dyn_alloc);
    return suite;
}

void RunAllTests(void) {
    
    printf("Running Test Suite: control\n");

    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * dyn = DynGetSuite();

    CuSuiteAddSuite(suite,dyn);
        
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);

    CuSuiteDelete(dyn);
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}
