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
    char fname[MAX_CHAR_LEN];
    
    strcpy(fname, "/Users/zrick/WORK/research_projects/GBM/triangulation");
    Triangulation tri(fname);   // call to libtriangulate

    tri.setAtlas(string("/Users/zrick/Work/research_projects/GBM/data/combined_atlas_labels.txt"));

    tri.writeGrid(3);
    std::cout << "Starting GBM \n";
    
    return 0;
}

