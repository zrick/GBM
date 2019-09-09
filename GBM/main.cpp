//
//  main.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "main.hpp"

int main(int argc, const char * argv[]) {
    int i;
    double dum=0.;
    string grid_file, grid_format, tri_file, atlas_file;
    
    Namelist nml("/Users/zrick/Work/research_projects/GBM/namelist.nml");

    tri_file  = nml.getVal_str("Files","triangulation");
    atlas_file= nml.getVal_str("Files","atlas");
    grid_file = nml.getVal_str("Files","grid");
    grid_format=nml.getVal_str("Files","grid_format");
    
    Triangulation tri((char *) &tri_file.c_str()[0]); // build triangulation
    tri.setAtlas(atlas_file);                         // add atlas to triangulation
    
    cout << "Total volume of Tetrahedrons: " << tri.volume  << '\n' ;
    cout << "Bounding Box:";
    for ( i=0; i<3; ++i)
        cout <<"[" << tri.box[2*i] << "," << tri.box[2*i+1] << "]" << ( i<2?" x ":"\n" );
    
    for ( i=0; i<tri.nVrt; ++i)
        if ( tri.vrt[i].vol > 0 )
            dum += tri.vrt[i].vol;
    cout <<"Total volume of Vertices: " << dum << '\n';

    tri.writeGrid(grid_file, grid_format);             // write grid to file

    
    std::cout << "Starting GBM \n";
    
    return 0;
}

