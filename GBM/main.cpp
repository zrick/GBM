//
//  main.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "main.hpp"
#include "triangulation_class.hpp"
int main(int argc, const char * argv[]) {
    // insert code here...
    char fname[MAX_CHAR_LEN];
    int status;
    
    
    strcpy(fname, "/Users/zrick/WORK/research_projects/GBM/triangles.txt");

    Triangulation tri(fname);
    
    std::cout << "Starting GBM \n";
 
    return 0;
}

