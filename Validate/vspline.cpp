//
//  vspline.cpp
//  vspline
//
//  Created by Cedrick Ansorge on 10.09.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "vspline.hpp"
int main(int argc, const char * argv[]) {
    
    double *dat_v, *dat_t;
    Triangulation *p_tri;

    Namelist nl("/Users/zrick/WORK/research_projects/GBM/vspline.nml");
    Triangulation tri((char *)&nl.getVal_str("Main","grid").c_str()[0]);
    p_tri=&tri;
    
    dat_v = (double *)malloc(tri.nVrt*sizeof(double));
    dat_t = (double *)malloc(tri.nTtr*sizeof(double));
    
    for (int i=0;i<p_tri->nVrt;++i)
        dat_v[i]=test_func(p_tri->vrt[i].p);
    
    exit(EXIT_SUCCESS);
    
}

double test_func(double *p){
    return (sin(p[0])*cos(p[1])*exp(p[2]));
}
