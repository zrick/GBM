//
//  main.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "main.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    int i;
    double dum=0.;
    char fname[MAX_CHAR_LEN];
    Namelist nml("/Users/zrick/Work/research_projects/GBM/namelist.nml");
    
    strcpy(fname, "/Users/zrick/WORK/research_projects/GBM/triangulation");
    Triangulation tri(fname);   // call to libtriangulate

    tri.setAtlas(string("/Users/zrick/Work/research_projects/GBM/data/atlas_all"));
    tri.writeGrid(3);
    
    cout << "Total volume of Tetrahedrons: " << tri.volume  << '\n' ;
    cout << "Bounding Box:";
    for ( i=0; i<3; ++i)
        cout <<"[" << tri.box[2*i] << "," << tri.box[2*i+1] << "]" << ( i<2?" x ":"\n" );
    
    for (i=0; i<tri.nVrt; ++i)
        dum += tri.vrt[i].vol;
    cout <<"Total volume of Vertices: " << dum << '\n';

    std::cout << "Starting GBM \n";
    
    return 0;
}

