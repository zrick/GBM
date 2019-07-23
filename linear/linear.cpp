//
//  linear.cpp
//  linear
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "linear.hpp"

using namespace std;

int qrdcmp(int nn)
{
    cout << "QR DECOMPOSITION NOT IMPLEMENTED\n";
    exit(EXIT_FAILURE);
    
    int sing=-1;
    return sing;
}

int qrslv(int nn)
{
    cout << "QR SOLVER NOT IMPLEMENTED\n";
    exit(EXIT_FAILURE);
    
    int status=-1;
    return status;
}

void crs3(double a[3], double b[3], double *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}

double dot3(double a[3], double b[3]){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void renorm3(double *v)
{
    double norm=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    for ( int i=0; i<3; ++i)
        v[i] /= norm;
    return;
}
