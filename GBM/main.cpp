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
    GBM gbm;
        
    vector<string> header;
    vector<string> labels; 
       

    //exit(EXIT_SUCCESS);
    
    GBMLog("====================================================\nSTARTING " + string(argv[0]) +" on " + gbmTime());
    
    if ( argc >1 ) {
        nl_name=string(argv[1]);
    }
    else   {
        nl_name=string("/Users/zrick/Work/research_projects/GBM/gbm.nml");
    }
    
    gbm.read_namelist(nl_name);
    gbm.init();
    gbm.tri.writeGrid(gbm.grid_file,gbm.grid_format);
  
    dum=0;
    for ( i=0; i<gbm.tri.nTtr; ++i)
        if ( gbm.tri.ttr[i].vol > 0 )
            dum += gbm.tri.ttr[i].vol;
    GBMLog("Total volume of Tetrahedrons: " +to_string(dum));
    
    dum=0.;
    for ( i=0; i<gbm.tri.nTri; ++i)
        dum+=gbm.tri.tri[i].area;
    GBMLog("Triangle Area: " +to_string(dum));
    
    dum=0;
    for ( std::vector<TriType *>::iterator tri_it = gbm.tri.hulTri.begin(); tri_it!= gbm.tri.hulTri.end(); ++tri_it) {
        dum+=(*tri_it)->area;
    }
    GBMLog("Surface Area: "+to_string(dum));
    
    gbm.tri.writeGrid("/Users/zrick/WORK/research_projects/GBM/test.vtu.xml","XML_VTK");
    gbm.tri.writePths("/Users/zrick/WORK/research_projects/GBM/pth.vtu.xml","XML_VTK");
    
    GBMLog("FINISHED GBM on " + gbmTime() + "====================================================");
      
    return 0;
}

