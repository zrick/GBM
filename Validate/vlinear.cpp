//
//  main.cpp
//  vlinear
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "vlinear.hpp"

int main(int argc, const char * argv[]) {
    
    int n,nm,n_block;
    double **a, **rhs,       ***am,***rhsm,            ***aO, ***rhsO;
    double *b, *c, *d, *lhs, **bm,**cm, **dm, **lhsm,  **bO, **cO, **dO, **lhsO;
    double dum,sum, max_res;
    int i,j,l,j_srt,j_end,sing;
    uint64_t t_srt;
    
    a=nullptr; b=nullptr; c=nullptr; d=nullptr;
    lhs=nullptr; rhs=nullptr;
    
    n=4;
    nm=571;
    n_block=int(n/2);
    
    cout << "ALLOCATING AND INITIALIZING ARRAYS\n";
    
    am  =(double***) malloc (n*sizeof(double **));
    rhsm=(double***) malloc (n*sizeof(double **));
    bm  =(double**)  malloc (n*sizeof(double *));
    cm  =(double**)  malloc (n*sizeof(double *));
    dm  =(double**)  malloc (n*sizeof(double *));
    lhsm=(double**)  malloc (n*sizeof(double *));
    
    for (i=0; i<n; ++i) {
        am[i]=(double**)malloc(n*sizeof(double*));
        rhsm[i]=(double**)malloc(n*sizeof(double*));
        
        bm[i]=(double *) malloc(nm*sizeof(double));
        cm[i]=(double *) malloc(nm*sizeof(double));
        dm[i]=(double *) malloc(nm*sizeof(double));
        lhsm[i]=(double*)malloc(nm*sizeof(double));
        
        for (l=0;l<nm; ++l){
            bm[i][l]=0.;
            cm[i][l]=0.;
            dm[i][l]=0.;
            lhsm[i][l]=0.;
        }
        
        for ( j=0; j<n; ++j) {
            am[i][j]  =(double *) malloc(nm*sizeof(double));
            rhsm[i][j]=(double *) malloc(nm*sizeof(double));
            for ( l=0; l<nm; ++l ) {
                am[i][j][l]=0.;
                rhsm[i][j][l]=0.;
            }
        }
    }
    
    aO  =(double***) malloc (nm*sizeof(double **));
    rhsO=(double***) malloc (nm*sizeof(double **));
    bO  =(double**)  malloc (nm*sizeof(double *));
    cO  =(double**)  malloc (nm*sizeof(double *));
    dO  =(double**)  malloc (nm*sizeof(double *));
    lhsO=(double**)  malloc (nm*sizeof(double *));
    
    for (l=0; l<nm; ++l) {
        aO[l]=(double**)malloc(n*sizeof(double*));
        rhsO[l]=(double**)malloc(n*sizeof(double*));
        
        bO[l]=(double *) malloc(n*sizeof(double));
        cO[l]=(double *) malloc(n*sizeof(double));
        dO[l]=(double *) malloc(n*sizeof(double));
        lhsO[l]=(double*)malloc(n*sizeof(double));
        
        for (i=0;i<n; ++i){
            bO[l][i]=0.;
            cO[l][i]=0.;
            dO[l][i]=0.;
            lhsO[l][i]=0.;
        }
        
        for ( i=0; i<n; ++i) {
            aO[l][i]  =(double *) malloc(nm*sizeof(double));
            rhsO[l][i]=(double *) malloc(nm*sizeof(double));
            for ( j=0; j<n; ++j ) {
                aO[l][i][j]=0.;
                rhsO[l][i][j]=0.;
            }
        }
    }

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
        for ( l=0; l<nm; ++l ){
            bm[i][l]=0.;
            bO[l][i]=0.;
            cm[i][l]=0.;
            cO[l][i]=0.;
   
        }
        
        for ( j=0; j<n; ++j) {
            a[i][j]=0.;
            rhs[i][j]=0.;
            for ( l=0; l<nm; ++l ) {
                am[i][j][l]=0;
                rhsm[i][j][l]=0.;
                
            }
        }
    }
    
    cout << "INITIALIZING \n";
    
    for (j=0; j<n; ++j) {
        b[j]=1.0 ;
        lhs[j]=1.0;
        a[j][j] =1.0;
        j_srt = j-n_block>0? j-n_block:0;
        j_end = j+n_block>n? n:j+n_block;
        
        for (l=0; l<nm; ++l) {
            bm[j][l]=1.0*sqrt(double(l));
            bO[l][j]=bm[j][l];
            lhsm[j][l]=1.0*sqrt(double(l));
            lhsO[l][j]=lhsm[j][l];
            am[j][j][l]=1.0*sqrt(double(l));
            aO[l][j][j]=am[j][j][l];
        }
        
        for (i=j_srt; i<j_end; ++i) {
            a[j][i]=1.-  fabs(double(i-j))/n_block;
            rhs[j][i]=1.-fabs(double(i-j))/n_block;
            for (l=0;l<nm;++l) {
                am[j][i][l] = a[j][i]*sqrt(l+1.);
                aO[l][j][i] = am[j][i][l];
                rhsm[j][i][l]=rhs[j][i]*sqrt(l+1.);
                rhsO[l][j][i]=rhsm[j][i][l];
            }
        }
    }
    
    cout <<"\nTesting QR decomposition for single-problem solver (n=" << n <<")\n";
    cout <<"============================================================\n";
    t_srt=current_time();
    sing=qrdcmp(a,n,n,c,d);
    cout << "qrdcmp time elapsed: "<<(current_time()-t_srt)/1000. << "ms \n";
    
    t_srt=current_time();
    qrsolv(a,n,n,c,d,b);
    cout <<"qrsolv time elapsed: "<<(current_time()-t_srt)/1000. << "ms \n";


    // Check single-problem solver
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
    cout <<"MAXIMUM RESIDUAL PER LINE (single-problem solver):" <<max_res<<'\n';
    
    // MULTI PROBLEM SOLVER WITH INNER LOOP
    cout <<"\nTesting QR decomposition for multi-problem solver; inner index (nm="<<nm<<")\n";
    cout <<"===========================================================================\n";
    t_srt=current_time();
    sing=m_in_qrdcmp(am,n,n,nm,cm,dm);
    cout << "m_in_qrdcmp time elapsed: "<<(current_time()-t_srt)/1000. << "ms \n";
    if ( sing > 0. ) {
        cout << "ERROR: Singularity encountered in qr decomposition \n";
        exit(EXIT_FAILURE);
    }
    
    t_srt=current_time();
    m_in_qrsolv(am,n,n,nm,cm,dm,bm);
    cout << "m_in_qrsolv time elapsed: "<<(current_time()-t_srt)/1000. << "ms \n";

    // Check multi-problem solver
    max_res=0.;
    for (l=0;l<nm;++l) {
        for(i=0;i<n;++i){
            sum=0.;
            for(j=0;j<n;++j)
                sum += rhsm[i][j][l]*bm[j][l];
            dum=fabs(sum-lhsm[i][l]);
            if ( dum > 1e-12 )
                cout <<"ERROR encountered in line " << i << "; Residual:" << dum << ' ' << lhsm[i][l] << '\n';
            max_res = dum > max_res? dum : max_res;
        }
    }
    cout <<"MAXIMUM RESIDUAL PER LINE (multi-problem solver):" << max_res <<"\n\n";


    // MULTI PROBLEM SOLVER WITH OUTER LOOP
    cout <<"\nTesting QR decomposition for multi-problem solver; outer index (nm="<<nm<<")\n";
    cout <<"===========================================================================\n";
    t_srt=current_time();
    sing=m_out_qrdcmp(aO,n,n,nm,cO,dO);
    cout << "m_out_qrdcmp time elapsed: "<<(current_time()-t_srt)/1000. << "ms \n";
    if ( sing > 0. ) {
        cout << "ERROR: Singularity encountered in qr decomposition \n";
        exit(EXIT_FAILURE);
    }
    
    t_srt=current_time();
    m_out_qrsolv(aO,n,n,nm,cO,dO,bO);
    cout << "m_out_qrsolv time elapsed: "<<(current_time()-t_srt)/1000. << "ms \n";
    
    // Check multi-problem solver
    max_res=0.;
    for (l=0;l<nm;++l) {
        for(i=0;i<n;++i){
            sum=0.;
            for(j=0;j<n;++j)
                sum += rhsO[l][i][j]*bO[l][j];
            dum=fabs(sum-lhsO[l][i]);
            if ( dum > 1e-12 )
                cout <<"ERROR encountered in line " << i << "; Residual:" << dum << ' ' << lhsO[l][i] << '\n';
            max_res = dum > max_res? dum : max_res;
        }
    }
    cout <<"MAXIMUM RESIDUAL PER LINE (multi-problem solver):" << max_res <<"\n\n";

    
    exit(EXIT_SUCCESS);
}
