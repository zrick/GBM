//
//  linear.hpp
//  linear
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef LINEAR_
#define LINEAR_

#include "constants.h"

int  qrdcmp(double **a,const int n,const int np,double *c,double *d);
void qrsolv(double **a, const int n, const int np, double *c, double *d, double *b);
void rsolv(double **a,const int n,const int np,double *d,double *b);

void crs3(double a[3], double b[3], double *c);
void renorm3(double *v);
double dot3(double a[3], double b[3]);

#endif
