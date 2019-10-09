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
    
    dum=0.;
    for ( i=0; i<gbm.tri.nTri; ++i)
        dum+=gbm.tri.tri[i].area;
     cout << "Triangle Area: " << dum;
    
    dum=0;
    for ( std::vector<TriType *>::iterator tri_it = gbm.tri.hul.begin(); tri_it!= gbm.tri.hul.end(); ++tri_it) {
        dum+=(*tri_it)->area;
    }
    cout << "; Surface Area: " << dum << std::endl;
    
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
    tri = Triangulation((char *) &(g.tri_file.c_str()[0] ) );
    if (g.use_atlas)
         tri.setAtlas(g.atlas_file);                         // add atlas to triangulation
     
    cout << "Total volume of Tetrahedrons: " << tri.volume  << '\n' ;
    cout << "Bounding Box:";
    for ( int i=0; i<3; ++i)
        cout <<"[" << tri.box[2*i] << "," << tri.box[2*i+1] << "]" << ( i<2?" x ":"\n" );
return;
}
