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

int qrdcmp(double **a,const int n,const int np,double *c,double *d){
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


int m_in_qrdcmp(double ***a,const int n,const int np,const int nm,double **c,double **d){
    int i,j,k,l,sgn[nm];
    double sum[nm],scale[nm],sigma[nm],tau[nm];
    int sing_all,sing[nm];
    
    sing_all=0;
    for (l=0;l<nm;++l)
        sing[l]=0;
        
    for (k=0; k<n-1; ++k) {
        for(l=0; l<nm; ++l)
            scale[l]=0.;
        for (i=k; i<n; ++i)
            for (l=0; l<nm; ++l){
                if ( a[k][i][l] > scale[l] || -a[k][i][l] > scale[l])
                    scale[l]= a[k][i][l] > 0? a[k][i][l] : -a[k][i][l];
            }
        for (l=0; l<nm; ++l)
            if ( scale[l] == 0. ) {
                sing[l] = 1;
                sing_all= 1;
                c[k][l]=0.;
                d[k][l]=0.;
        } else {
            
            for ( i=k; i<n; ++i)
                for (l=0;l<nm;++l)
                a[k][i][l]=a[k][i][l]/scale[l];
            
            for (l=0; l<nm; ++l)
                sum[l]=0. ;
            
            for ( i=k; i<n; ++i)
                for (l=0;l<nm;++l)
                    sum[l] = sum[l]+a[k][i][l]*a[k][i][l];
            
            for (l=0; l<nm; ++l){
                sgn[l] = a[k][k][l] >=0. ? 1 : -1;
                sigma[l]=sgn[l] * sqrt(sum[l]);
                a[k][k][l]=a[k][k][l]+sigma[l] ;
                c[k][l]=sigma[l]*a[k][k][l];
                d[k][l]=-scale[l]*sigma[l];
            }
            
            for ( j=k; j<n; ++j ){
                for (l=0; l<nm; ++l)
                    sum[l]=0;
                for (i=k; i<n; ++i )
                    for (l=0; l<nm; ++l)
                        sum[l]=sum[l]+a[k][i][l]*a[j][i][l];
                for (l=0; l<nm; ++l)
                    tau[l]=sum[l]/c[k][l] ;
                for (i=k; i<n; ++i)
                    for (l=0; l<nm; ++l)
                        a[j][i][l]=a[j][i][l]-tau[l]*a[k][i][l];
            }
        }
    }

    for (l=0; l<nm; ++l) {
        d[n-1][l]=a[n-1][n-1][l];
        if ( d[n-1][l] == 0.)
            sing_all=1;
    }
    return sing_all;
}

void m_in_qrsolv(double ***a, const int n, const int np, const int nm, double **c, double **d, double **b){
    int i,j,l;
    double sum[nm],tau[nm];
    
    for (j=0; j<n-1; ++j ) {
        for (l=0;l<nm;++l)
            sum[l]=0. ;
        for (i=j; i<n; ++i)
            for (l=0;l<nm;++l)
                sum[l]+=a[j][i][l]*b[i][l];
        for (l=0;l<nm;++l)
            tau[l]=sum[l]/c[j][l] ;
        for (i=j; i<n; ++i)
            for (l=0;l<nm;++l)
                b[i][l]=b[i][l]-tau[l]*a[j][i][l];
    }
    m_in_rsolv(a,n,np,nm,d,b) ;
    
    return;
}

void m_in_rsolv(double ***a,const int n,const int np, const int nm, double **d,double **b)
{
    int i,j,l;
    double sum[nm];
    
    for (l=0;l<nm;++l)
        b[n-1][l]=b[n-1][l]/d[n-1][l] ;
    for (i=n-2; i>=0; --i) {
        for (l=0;l<nm;++l)
            sum[l]=0;
        for (j=i+1; j<n; ++j)
            for (l=0;l<nm;++l)
                sum[l]+= a[j][i][l]*b[j][l] ;
        for (l=0;l<nm;++l)
            b[i][l] = (b[i][l]-sum[l])/d[i][l];
    }
    
    return;
}

int m_out_qrdcmp(double ***a,const int n,const int np,const int nm, double **c,double **d){
    int i,j,k,l,sgn;
    double sum,scale,sigma,tau;
    int sing=0 ;
    
    for (l=0; l<nm; ++l){
    for (k=0; k<n-1; ++k) {
        scale=0.;
        for (i=k; i<n; ++i)
            if ( a[l][k][i] > scale || -a[l][k][i] > scale)
                scale= a[l][k][i] > 0? a[l][k][i] : -a[l][k][i];
        
        if ( scale == 0. ) {
            sing = 1 ;
            c[l][k]=0. ;
            d[l][k]=0. ;
        } else {
            
            for ( i=k; i<n; ++i)
                a[l][k][i]=a[l][k][i]/scale;
            
            sum=0. ;
            
            for ( i=k; i<n; ++i)
                sum = sum+a[l][k][i]*a[l][k][i];
            
            sgn = a[l][k][k] >=0. ? 1 : -1;
            sigma=sgn * sqrt(sum);
            a[l][k][k]=a[l][k][k]+sigma ;
            c[l][k]=sigma*a[l][k][k];
            d[l][k]=-scale*sigma;
            
            for ( j=k; j<n; ++j ){
                sum=0;
                for (i=k; i<n; ++i )
                    sum=sum+a[l][k][i]*a[l][j][i];
                tau=sum/c[l][k] ;
                for (i=k; i<n; ++i)
                    a[l][j][i]=a[l][j][i]-tau*a[l][k][i];
            }
        }
    }
    
    d[l][n-1]=a[l][n-1][n-1];
    if ( d[l][n-1] == 0.)
        sing=1 ;
    }
    
    return sing;
}

void m_out_qrsolv(double ***a, const int n, const int np, const int nm, double **c, double **d, double **b)
{
    int i,j,l;
    double sum,tau;
    
    for (l=0; l<nm; ++l) {
    for (j=0; j<n-1; ++j ) {
        sum=0. ;
        for (i=j; i<n; ++i)
            sum+=a[l][j][i]*b[l][i];
        tau=sum/c[l][j] ;
        for (i=j; i<n; ++i)
            b[l][i]=b[l][i]-tau*a[l][j][i];
    }}
    m_out_rsolv(a,n,np,nm,d,b) ;
    
    return;
}

void m_out_rsolv(double ***a,const int n,const int np, const int nm, double **d,double **b) {
    int i,j,l;
    double sum;
    
    for (l=0;l<nm;++l) {
    b[l][n-1]=b[l][n-1]/d[l][n-1] ;
    for (i=n-2; i>=0; --i) {
        sum=0;
        for (j=i+1; j<n; ++j)
            sum+= a[l][j][i]*b[l][j] ;
        b[l][i] = (b[l][i]-sum)/d[l][i];
    }}
    
    return;
}



double dist3(double a[3], double b[3])
{
    double d=0.;
    for (int i=0; i<3; ++i )
        d += (b[i]-a[i])*(b[i]-a[i]);
    return(d);
}

void crs3(double a[3], double b[3], double *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return;
}

double dot3(double a[3], double b[3]){
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

void renorm3(double *v)
{
    double norm=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    for ( int i=0; i<3; ++i)
        v[i] /= norm;
    return;
}

double tetraVolume(double a[3], double b[3], double c[3], double d[3]) {
// dot ( cross(b-a,c-a), d-a ) / 6
//                   (b-a)            x   (c-a)                     (d-a)
   return fabs( ( (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]) ) * (d[0]-a[0])
           +( (b[2]-a[2])*(c[0]-a[0]) - (b[2]-a[0])*(c[1]-a[2]) ) * (d[1]-a[1])
           +( (b[0]-a[0])*(c[1]-a[1]) - (b[2]-a[1])*(c[1]-a[0]) ) * (d[2]-a[2]) ) /6.;
}
