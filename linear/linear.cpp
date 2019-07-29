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

int qrdcmp(double **a,int n,int np,double *c,double *d){
    int i,j,k,sgn;
    double sum,scale,sigma,tau;
    int sing=0 ;
    
    for (k=0; k<n-1; ++k) {
        scale=0.;
        for (i=k; i<n; ++i)
            if ( a[k][i] > scale || -a[k][i] > scale)
                scale= a[k][i] > 0? a[k][i] : -a[k][i];
        
        if ( scale == 0. ) {
            sing = 1 ;
            c[k]=0. ;
            d[k]=0. ;
        } else {
            
            for ( i=k; i<n; ++i)
                a[k][i]=a[k][i]/scale;
            
            sum=0. ;
            
            for ( i=k; i<n; ++i)
                sum = sum+a[k][i]*a[k][i];
            
            sgn = a[k][k] >=0. ? 1 : -1;
            sigma=sgn * sqrt(sum);
            a[k][k]=a[k][k]+sigma ;
            c[k]=sigma*a[k][k];
            d[k]=-scale*sigma;
            
            for ( j=k; j<n; ++j ){
                sum=0;
                for (i=k; i<n; ++i )
                    sum=sum+a[k][i]*a[j][i];
                tau=sum/c[k] ;
                for (i=k; i<n; ++i)
                    a[j][i]=a[j][i]-tau*a[k][i];
            }
        }
    }
    
    d[n-1]=a[n-1][n-1];
    if ( d[n-1] == 0.)
        sing=1 ;
    
    return sing;
}

void qrsolv(double **a, const int n, const int np, double *c, double *d, double *b)
{
    int i,j;
    double sum,tau;
    
    for (j=0; j<n-1; ++j ) {
        sum=0. ;
        for (i=j; i<n; ++i)
            sum+=a[j][i]*b[i];
        tau=sum/c[j] ;
        for (i=j; i<n; ++i)
            b[i]=b[i]-tau*a[j][i];
    }
    rsolv(a,n,np,d,b) ;
    
    return;
}

void rsolv(double **a,const int n,const int np,double *d,double *b)
{
    int i,j;
    double sum;
    
    b[n-1]=b[n-1]/d[n-1] ;
    for (i=n-2; i>=0; --i) {
        sum=0;
        for (j=i+1; j<n; ++j)
            sum+= a[j][i]*b[j] ;
        b[i] = (b[i]-sum)/d[i];
    }
    
    return;
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
