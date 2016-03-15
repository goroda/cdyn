// Copyright (c) 2016, Massachusetts Institute of Technology
//
// Author: Alex A. Gorodetsky 
// Contact: goroda@mit.edu

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <stdlib.h>

struct Drift;
struct Drift * 
drift_alloc_controlled(size_t, size_t,
                       int (*)(double,const double*,const double*,
                               double*,double*,void*),
                       void *);
struct Drift * 
drift_alloc_uncontrolled(size_t,
                         int (*)(double,const double*,double*,void*),
                         void *);                                                
void drift_free(struct Drift *);
int drift_eval(const struct Drift *,double,const double *,
               const double *, double *, double *);

////////////////////////////////////
struct Diff;
struct Diff * 
diff_alloc_controlled(size_t,size_t,size_t,
                      int (*)(double,const double*,const double*,
                              double*,double*,void*),
                      void *);
struct Diff * 
diff_alloc_uncontrolled(size_t, size_t,
                        int (*)(double,const double*,double*,void*),
                        void *);
void diff_free(struct Diff *);
int diff_eval(const struct Diff *,double,const double*,const double*,
              double*,double*);


///////////////////////////////////////////////////////////

struct Dyn * dyn_alloc(struct Drift *, struct Diff *);
void dyn_free(struct Dyn *);
void dyn_free_deep(struct Dyn *);

size_t dyn_get_dx(const struct Dyn *);
size_t dyn_get_dw(const struct Dyn *);
size_t dyn_get_du(const struct Dyn *);
int dyn_eval(const struct Dyn *,double,const double*,const double*,double*,
             double*,double*,double*);


////////////////////////////////////////////////////////////

#endif
