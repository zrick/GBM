//
//  main.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "main.hpp"
#include "gbm_data.h"

int main(int argc, const char * argv[]) {
    int i;
    double dum;
    string nl_name;
    GBM_Data gbm;

    
    if ( argc >1 ) {
        nl_name=string(argv[1]);
    }
    else   {
        nl_name=string("/Users/zrick/Work/research_projects/GBM/namelist.nml");
    }
    
    gbm_read_namelist(nl_name,gbm);
    gbm_init(gbm,gbm.tri);
        
     
    std::cout << "Starting GBM \n";

    gbm.tri.writeGrid(gbm.grid_file,gbm.grid_format);
  
    dum=0;
    for ( i=0; i<gbm.tri.nVrt; ++i)
        if ( gbm.tri.vrt[i].vol > 0 )
            dum += gbm.tri.vrt[i].vol;
    cout <<"Total volume of Vertices: " << dum << '\n';
    
    std::cout <<"GBM FINISHED" << std::endl;
    
    return 0;
}

void gbm_read_namelist(string &nl_file,GBM_Data &g){
    Namelist *nl;
    
    g.nml=Namelist(nl_file);
    nl=&(g.nml);
    
    // Mandatory arguments -- no default if missing from namelist.
    
    g.tri_file  = nl->getVal_str("Files","triangulation");
    g.grid_file = nl->getVal_str("Files","grid");
    g.grid_format=nl->getVal_str("Files","grid_format");
    
    g.use_atlas=false;
    if ( nl->hasVal("Files","atlas") )
    {
        g.use_atlas=true;
        g.atlas_file= nl->getVal_str("Files","atlas");
    }
    
    return;
}

void gbm_init(GBM_Data &g, Triangulation &tri){
    double dum=0;
    double *p1, *p2, *p3, *p4;
    double d1[3], d2[3], d3[3], d4[3], d5[3];
    int i;

    tri = Triangulation((char *) &(g.tri_file.c_str()[0] ) );
    if (g.use_atlas)
         tri.setAtlas(g.atlas_file);                         // add atlas to triangulation
     
    for (i=0; i<tri.nTtr; ++i)
        if ( 6.*4.*tri.ttr[i].vol != 1 )
        {
            cout << tri.ttr[i].vol;
            p1 = & (tri.vrt[tri.ttr[i].vrt[0]].p[0]);
            p2 = & (tri.vrt[tri.ttr[i].vrt[1]].p[0]);
            p3 = & (tri.vrt[tri.ttr[i].vrt[2]].p[0]);
            p4 = & (tri.vrt[tri.ttr[i].vrt[3]].p[0]);
            
            for (int id=0; id<3; ++id){
                d1[id]=p2[id]-p1[id];
                d2[id]=p3[id]-p1[id];
                d3[id]=p4[id]-p1[id];
            }
            crs3(d1,d2,&d4[0]);
            /*cout << "(" << p1[0] <<","<<p1[1]<<","<<p1[2]<<"),";
            cout << "(" << p2[0] <<","<<p2[1]<<","<<p2[2]<<"),";
            cout << "(" << p3[0] <<","<<p3[1]<<","<<p3[2]<<"),";
            cout << "(" << p4[0] <<","<<p4[1]<<","<<p4[2]<<")" << std::endl << "   ";
            cout << "(" << d1[0] <<","<<d1[1]<<","<<d1[2]<<"),";
            cout << "(" << d2[0] <<","<<d2[1]<<","<<d2[2]<<"),";
            cout << "(" << d3[0] <<","<<d3[1]<<","<<d3[2]<<"),";
            cout << "(" << d4[0] <<","<<d4[1]<<","<<d4[2]<<"),";*/
           
            
            
            cout <<" " << fabs(dot3(d4,d3)/6.) <<  std::endl;
        }
    cout << std::endl;
    cout << "Total volume of Tetrahedrons: " << tri.volume  << '\n' ;
    cout << "Bounding Box:";
    for ( i=0; i<3; ++i)
        cout <<"[" << tri.box[2*i] << "," << tri.box[2*i+1] << "]" << ( i<2?" x ":"\n" );
     
    for ( i=0; i<tri.nVrt; ++i)
        if ( tri.vrt[i].vol > 0 )
            dum += tri.vrt[i].vol;
    cout <<"Total volume of Vertices: " << dum << '\n';
return;
}
