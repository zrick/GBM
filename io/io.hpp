//
//  io.hpp
//  io
//
//  Created by Cedrick Ansorge on 30.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef io_
#define io_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include "constants.h"
#include "utils.hpp"
#include "types.h"

// #include "namelist.hpp"

/* The classes below are exported */
#pragma GCC visibility push(default)

using namespace std;

template<typename T> // explicit instantiation to types <double, int> in io.cpp
void vtkXMLWriteDataArray(ofstream &gfile, vector<string> att, int nval, T *val);

void vtkXMLFileOpen(string &name,int np,int nc,ofstream &f);
void vtkXMLFileClose(ofstream &f);

#pragma GCC visibility pop
#endif
