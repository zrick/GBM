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
    Triangulation(char *fname,bool p[3]);
 
    // Function Members
    void writeGrid(string grid_file, string grid_format);
    void printHull(int format);
    
    void printVertex(VertexType *v, int level);
    void printEdge(EdgeType *e, int level);
    void printTriangle(TriType *t, int level);
    void printTetrahedron(TetraType *t, int level);
    
    bool hasAtlas();
    void setAtlas(string fname);
    
    void ConstructHull();
    void ConstructHalo();

    void TtrSplines_lhs();
    void TtrSplines_rhs(double *data);
    void CentroidSplines(double *v); 
    
    // Data Members
    int nDim, nVrt, nEdg, nTri, nTtr, nHul;
    int nTtr_inner, nVrt_inner;
    double box[6], halo_box[6];
    double volume;
    
    vector<TetraType>  ttr; // 3D Tetrahedrons
    vector<VertexType> vrt; // 0D Vertex Points
    vector<TriType>    tri; // 2D Triangles (surfaces of Tetrahedrons)
    vector<EdgeType>   edg; // 1D Edges of Tetrahedrons and Triangles
    //  
    vector<TriType *>  hul; // Exterior Surface of Triangulation
    vector<string>   atlas; // Names of atlas regions; first element is atlas file
    
private:
    void read_mesh(char *fname);
    void construct_normal(TetraType *t, TriType tri_loc, int ittr_loc);

    int addVertex(istringstream &l,bool halo=false);
    int addVertex(double p[3],     bool halo=false);
    int addTetra(istringstream &l, bool halo=false);
    int addTetra(int p[4],         bool halo=false);
    int addEdges(int it, TetraType *t);
    int addTris(int it, TetraType *t);
    void addHalo();

    int isTetra(int p[4], int start=0, int end=-1);
    int isTriangle(int v[3]);
    int isEdge(int v0, int v1);
    int isVertex(double p[3], int start=0, int end=-1);
    
    void ttrPeriodic(int it, int iv, int idim, int *ttr_new);
    void ttrPeriodic(int it, int iv2, bool per[3], int *ttr_new);

    double ***QRSplines_a;
    double **QRSplines_b;
    double **QRSplines_c;
    double **QRSplines_d;
    double *QRSplines_bare;
    
    bool periodic[3]; 
};

#endif // ifndef triangulate_h
