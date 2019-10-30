//
//  main.hpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef main_h

#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <utility>
#include <algorithm>

#include "utils.hpp"
#include "triangulate.hpp"
#include "constants.h"
#include "io.hpp"
#include "types.h"

class GBM{
    
public:
    Namelist nml;
    Triangulation tri;

    string grid_file;
    string grid_format;
    string atlas_file;
    string pathes_file;
    string tau_file; 
    string tri_file;

    bool use_atlas;
    bool use_pathes;
    bool tau_from_file;
    bool periodic[3];

    double domain[6];
    
    double *tauVrt;
    
    void read_namelist(string &nl_file);
    void init();
    
    
private:
    friend Triangulation;
    int tauLine; 
};

#define main_h
#endif /* main_h */
