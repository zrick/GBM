//
//  main.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "main.hpp"
#include "utils.hpp"
#include "triangulate.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    char fname[MAX_CHAR_LEN];
    int a[4]={0,3,2,1};
    strcpy(fname, "/Users/zrick/WORK/research_projects/GBM/triangulation");

    
    sort4(a); // test call to libutils.
    
    Triangulation tri(fname);   // call to libtriangulate
    //Triangulation tri;
    
    tri.printGrid(0);
    
    std::cout << "Starting GBM \n";
    
    return 0;
}

