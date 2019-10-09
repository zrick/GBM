//
//  vspline.hpp
//  vspline
//
//  Created by Cedrick Ansorge on 10.09.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef vspline_hpp
#define vspline_hpp

#include <iostream>
#include <stdio.h>
#include <cstring> 

#include "io.hpp"
#include "utils.hpp"
#include "triangulate.hpp"
#include "linear.hpp"

using namespace std;

double test_func(double *p, double *k);
void writeSpline_test(string grid_file, Triangulation *tri, double *dat_v, double *dat_t);

#endif /* vspline_hpp */
