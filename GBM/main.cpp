//
//  main.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "triangulation_class.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    char fname[MAX_CHAR_LEN];
    
    strcpy(fname, "/Users/zrick/WORK/research_projects/GBM/triangulation");

    Triangulation tri(fname);
    //Triangulation tri;
    
    tri.printGrid(0);
    
    std::cout << "Starting GBM \n";
    
    return 0;
}

