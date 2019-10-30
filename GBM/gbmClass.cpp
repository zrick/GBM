//
//  gbmClass.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 30.10.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "main.hpp"

void GBM::read_namelist(string &nl_file){
    Namelist *nl;
    
    nml=Namelist(nl_file);
    nl=&(nml);
    
    // mandatory arguments -- no default if missing from namelist.
    
    tri_file  = nl->getVal_str("Files","triangulation");
    grid_file = nl->getVal_str("Files","grid");
    grid_format=nl->getVal_str("Files","grid_format");
    
    // Optional values -- defaults set.
    
    use_atlas=false;
    if ( nl->hasVal("Files","atlas") ) {
        use_atlas=true;
        atlas_file= nl->getVal_str("Files","atlas");
    }
    
    use_pathes=false;
    if ( nl->hasVal("Files","pathes") ) {
        use_pathes=true;
        pathes_file=nl->getVal_str("Files","pathes");
    }
    
    tau_from_file=false;
    if ( nl->hasVal("Files","tauFile") ) {
        tau_from_file=true;
        tau_file=nl->getVal_str("Files","tauFile");
        tauLine=1;
        if ( nl->hasVal("Files","tauLine") )
            tauLine=nl->getVal_int("Files","tauLine");
    }
    
    for(int i=0; i<3; ++i)   // Default is non-periodic
        periodic[i]=false;
    if ( nl->hasVal("Grid","periodic_x" ) )
        periodic[0]=nl->getVal_bool("Grid","periodic_x");
    if ( nl->hasVal("Grid","periodic_y" ) )
        periodic[1]=nl->getVal_bool("Grid","periodic_y");
    if ( nl->hasVal("Grid","periodic_z" ) )
        periodic[2]=nl->getVal_bool("Grid","periodic_z");
    
    nl->getList_dbl("Grid", "domain", domain, 6);
    
    return;
}

void GBM::init(){
    // initialize triangulation
    tri = Triangulation((char *) &(tri_file.c_str()[0] ),periodic,domain);
    string str;
    stringstream sstr;

    GBMLog("Total volume of Tetrahedrons: " +to_string(tri.volume));
    sstr << "Bounding Box:" << setprecision(4) << scientific;
    for ( int i=0; i<3; ++i)
        sstr <<  "[" << tri.box[2*i] << "," << tri.box[2*i+1] << "]" << ( i<2?" x ": "" );
    GBMLog(sstr.str());
    
    if (use_atlas)
        tri.setAtlas(atlas_file);

    if ( periodic[0] or periodic[1] or periodic[2] )
        tri.ConstructHalo();

    if ( use_pathes )
        tri.ConstructPathes(pathes_file);
    
    if ( tau_from_file ) {
        string dum;
        tauVrt = (double *) malloc ( tri.nVrt *sizeof(double) );
        readCSVLine(tau_file, tauLine, tauVrt, dum);
    }
    
    return;
}
