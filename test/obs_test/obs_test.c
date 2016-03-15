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
#include "observations.h"

size_t dx = 3;
size_t dy = 2;

int obsf(double time, const double * x, double * out, void * arg)
{
    double coeff = *(double *)arg;
    out[0] = exp(-time *10)*x[1]*x[2];
    out[1] = coeff + sin(x[0] + x[2]);
    return 0;
}

////////////////////////////////////////////////////////////////////

void Test_observation_alloc(CuTest * tc)
{

    fprintf(stdout,"Testing allocation and evaluation of Observation\n");
    double coeff = 2.0;
    struct Observation * obs = observation_alloc(dx,dy,obsf,&coeff);

    double test_pt[3] = {0.3, 2.0, -0.4};
    double time = 2.3;
    double obsvals[2];
    double obsvals2[2];

    int res = observation_observe(obs,time,test_pt,obsvals);
    CuAssertIntEquals(tc,0,res);

    obsf(time,test_pt,obsvals2,&coeff);
    for (size_t ii = 0; ii < dy; ii++){
        CuAssertDblEquals(tc,obsvals2[ii],obsvals[ii],1e-20);
    }

    observation_free(obs); obs = NULL;
}

CuSuite * ObsGetSuite(){

    CuSuite * suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, Test_observation_alloc);
    return suite;
}

void RunAllTests(void) {
    
    printf("Running Test Suite: control\n");

    CuString * output = CuStringNew();
    CuSuite * suite = CuSuiteNew();
    
    CuSuite * obss = ObsGetSuite();

    CuSuiteAddSuite(suite,obss);
        
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s \n", output->buffer);

    CuSuiteDelete(obss);
    CuStringDelete(output);
    free(suite);
   
}

int main(void) {
    RunAllTests();
}
