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

// #include "namelist.hpp"

/* The classes below are exported */
#pragma GCC visibility push(default)

using namespace std;
class Attribute{
public:
    // Constructors
    Attribute(string &name, string &val);
    string attName,attVal;
};



class Group{
public:
    // Constructors
    Group(string &s);
    string name;
    void addAttribute(string n, string v);
    void getAttribute(string &n, string &v);
    
    std::vector<Attribute> att;
private:
    // Functions
    
    // Members

};


class Namelist{
    
public:
    // Constructors
    Namelist();
    Namelist(char *fname);
    Namelist(const char fname[]);
    
    int    getVal_int(string group,string name);
    double getVal_dbl(string group,string name);
    string getVal_str(string group,string name);
    bool   getVal_bool(string group,string name);
    
private:
    //Functions
    void read_namelist(string s);
    string getValStr(string group,string var);

    // Members
    std::vector<Group> grp;
    std::vector<string>grpNames;
    
    
};

template<typename T> // explicit instantiation to types <double, int> in io.cpp
void vtkXMLWriteDataArray(ofstream &gfile, vector<string> att, int nval, T *val);

void vtkXMLFileOpen(string &name,int np,int nc,ofstream &f);
void vtkXMLFileClose(ofstream &f);

#pragma GCC visibility pop
#endif
