//
//  main.cpp
//  vlinear
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "vlinear.hpp"

int main(int argc, const char * argv[]) {
    
    int n,n_block;
    double **a, **rhs;
    double *b, *c, *d, *lhs;
    double dum,sum, max_res;
    int i,j,j_srt,j_end,sing;
    uint64_t t_srt;
    
    a=nullptr; b=nullptr; c=nullptr; d=nullptr;
    lhs=nullptr; rhs=nullptr;
    
    n=100;
    n_block=int(n/2);
    
    cout << "ALLOCATING AND INITIALIZING ARRAYS\n";
    a  =(double**) malloc (n*sizeof(double *));
    rhs=(double**) malloc (n*sizeof(double *));
    b  =(double*)  malloc (n*sizeof(double));
    c  =(double*)  malloc (n*sizeof(double));
    d  =(double*)  malloc (n*sizeof(double));
    lhs=(double*)  malloc (n*sizeof(double));
    for (i=0; i<n; ++i) {
        a[i]=(double*)malloc(n*sizeof(double));
        rhs[i]=(double*)malloc(n*sizeof(double));
        
        b[i]=0.;
        c[i]=0.;
        d[i]=0.;
        lhs[i]=0.;
        for ( j=0; j<n; ++j) {
            a[i][j]=0;
            rhs[i][j]=0.;
        }
    }
    
    cout << "INITIALIZING \n";
    
    for (j=0; j<n; ++j) {
        b[j]=1.0 ;
        lhs[j]=1.0;
        a[j][j] =1.0;
        j_srt = j-n_block>0? j-n_block:0;
        j_end = j+n_block>n? n:j+n_block;

        for (i=j_srt; i<j_end; ++i) {
            a[j][i]=1.-  fabs(double(i-j))/n_block;
            rhs[j][i]=1.-fabs(double(i-j))/n_block;
        }
    }
    
    cout <<"Testing QR decomposition\n";
    t_srt=current_time();
    sing=qrdcmp(a,n,n,c,d);
    cout << "qrdcmp time elapsed: "<<current_time()-t_srt << "micro s \n";
    
    if ( sing > 0. ) {
        cout << "ERROR: Singularity encountered in qr decomposition \n";
        exit(EXIT_FAILURE);
    }

    t_srt=current_time();
    qrsolv(a,n,n,c,d,b);
    cout <<"qrsolv time elapsed: "<<current_time()-t_srt << "micro s \n";
    
    max_res=0.;
    for(i=0;i<n;++i)
    {
        sum=0.;
        for(j=0;j<n;++j)
            sum += rhs[i][j]*b[j];
        dum=fabs(sum-lhs[i]);
        if ( dum > 1e-12 )
            cout <<"ERROR encountered in line " << i << "; Residual:" << dum << ' ' << lhs[i] << '\n';
        max_res = dum > max_res? dum : max_res;
    }

    cout <<"MAXIMUM RESIDUAL PER LINE:" <<max_res<<'\n';
 
    /*for (j=0;  j<n; ++j) {
        for ( i=0; i<n; ++i ) {
            cout << rhs[j][i] <<" " ;
        }
        cout << '=' << lhs[j] << '\n';
    }
    cout <<"Solution: ";
    for (j=0; j<n; ++j) {
        cout << b[j] << " ";
    }
    cout << '\n';
     */
    exit(EXIT_SUCCESS);
}
