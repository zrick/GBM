//
//  triangulate.hpp
//  triangulate
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef triangulate_h
#define triangulate_h

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <utility>
#include <vector>

#include "constants.h"
#include "types.h"

#include "linear.hpp"
#include "utils.hpp"
#include "io.hpp"


using namespace std;

class Triangulation{
    
public:
    // Constructors
    Triangulation();
    Triangulation(char *fname);
 
    // Function Members
    void writeGrid(string grid_file, string grid_format);
    void printHull(int format);
    
    void printVertex(VertexType *v, int level);
    void printEdge(EdgeType *e, int level);
    void printTriangle(TriType *t, int level);
    void printTetrahedron(TetraType *t, int level);
    bool hasAtlas();
    void setAtlas(string fname);
    
    // Data Members
    int nDim, nVrt, nEdg, nTri, nTtr, nHul;
    double box[6]; 
    double volume;
    
    TetraType *ttr;              // 3D Tetrahedrons
    VertexType *vrt;             // 0D Vertex Points
    std::vector<TriType>    tri; // 2D Triangles (surfaces of Tetrahedrons)
    std::vector<EdgeType>   edg; // 1D Edges of Tetrahedrons and Triangles
    std::vector<TriType *>  hul; // Exterior Surface of Triangulation

    std::vector<string>   atlas; // Names of atlas regions; first element is atlas file
    
private:
    void read_mesh(char *fname);
    void construct_normal(TetraType *t, TriType tri_loc, int ittr_loc);

    int add_edges(int it, TetraType *t);
    int add_tris(int it, TetraType *t);
    int is_triangle(int v[3]);
    int is_edge(int v0, int v1);
};

#endif // ifndef triangulate_h
