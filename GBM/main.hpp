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

typedef struct GBM_Data{
    string grid_file;
    string grid_format;
    string atlas_file;
    bool use_atlas; 
    string tri_file;
    Namelist nml;
    Triangulation tri;
    bool periodic[3];
    double domain[6];
} GBM_Data;

void gbm_read_namelist(string &nl_file, GBM_Data &g);
void gbm_init(GBM_Data &g, Triangulation &t);
#define main_h
#endif /* main_h */
