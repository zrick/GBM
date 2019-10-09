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
    double dum=0.;
    string nl_name;
    GBM_Data gbm;
    GBM_Data *p_gbm;
    
    p_gbm=&gbm;

    if ( argc >1 ) {
        nl_name=string(argv[1]);
    }
    else   {
        nl_name=string("/Users/zrick/Work/research_projects/GBM/namelist.nml");
    }
    
    gbm_read_namelist(nl_name,p_gbm);
    gbm_init(p_gbm);
    
    cout << "Total volume of Tetrahedrons: " << p_gbm->tri->volume  << '\n' ;
    cout << "Bounding Box:";
    for ( i=0; i<3; ++i)
        cout <<"[" << p_gbm->tri->box[2*i] << "," << p_gbm->tri->box[2*i+1] << "]" << ( i<2?" x ":"\n" );
    
    dum=0;
    for ( i=0; i<p_gbm->tri->nVrt; ++i)
        if ( p_gbm->tri->vrt[i].vol > 0 )
            dum += p_gbm->tri->vrt[i].vol;
    cout <<"Total volume of Vertices: " << dum << '\n';
    
    std::cout << "Starting GBM \n";

    p_gbm->tri->writeGrid(p_gbm->grid_file,p_gbm->grid_format);
    
    return 0;
}

void gbm_read_namelist(string &nl_file,GBM_Data *g){
    Namelist *nl;
    
    g->nml=Namelist(nl_file);
    nl=&g->nml;
    
    // Mandatory arguments -- no default if missing from namelist.
    
    g->tri_file  = nl->getVal_str("Files","triangulation");
    g->grid_file = nl->getVal_str("Files","grid");
    g->grid_format=nl->getVal_str("Files","grid_format");
    
    g->use_atlas=false;
    if ( nl->hasVal("Files","atlas") )
    {
        g->use_atlas=true;
        g->atlas_file= nl->getVal_str("Files","atlas");
    }

    
    
    return;
}

void gbm_init(GBM_Data *g){
    
    g->tri=(Triangulation *) malloc(sizeof(Triangulation));
    g->tri[0]= Triangulation((char *) &(g->tri_file.c_str()[0]));// build triangulation

    if (g->use_atlas)
        g->tri->setAtlas(g->atlas_file);                         // add atlas to triangulation
    
    g->tri->writeGrid(g->grid_file, g->grid_format);             // write grid to file

    return;
}
