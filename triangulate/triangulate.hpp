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
#include <iterator>
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
    Triangulation(char *fname,bool p[3], double *dom);
 
    // Function Members
    int allocateWrk(); 
    void writeGrid(string grid_file, string grid_format);
    void writePths(string path_file, string grid_format);
    void printHull(int format);
    
    void printVertex(VertexType *v, int level);
    void printEdge(EdgeType *e, int level);
    void printTriangle(TriType *t, int level);
    void printTetrahedron(TetraType *t, int level);
    
    bool hasAtlas();
    bool hasHalo();
    void setAtlas(string fname);
    
    void ConstructHull();
    void ConstructHalo();
    void ConstructPathes(string pfile); 
    
    void TtrSplinesAlloc(int sType);
    void TtrSplinesLHS(int sType);
    void TtrSplinesRHS(int sType, double *data);
    void TtrCentroidSplines(int sType, double *v); 
    void TtrDerivativeSplines(int sType, int direction, double *v);

    // Data Members
    int nDim, nVrt, nEdg, nTri, nTtr, nHul, nPth;
    int nTtr_inner, nVrt_inner;
    double box[6], halo_box[6];
    double volume;
    
    vector<TetraType>  ttr; // 3D Tetrahedrons
    vector<VertexType> vrt; // 0D Vertex Points
    vector<TriType>    tri; // 2D Triangles (surfaces of Tetrahedrons)
    vector<EdgeType>   edg; // 1D Edges of Tetrahedrons and Triangles
    //  
    vector<TriType *>    hulTri; // Exterior Surface Triangles of Triangulation
    vector<int>          hulVrt; // Exterior Surface Vertices of Triangulation 
    vector<string>   atlas; // Names of atlas regions; first element is atlas file
    
    HashBox hbox[GBM_NBOX+2][GBM_NBOX+2][GBM_NBOX+2];
    double hash_dx[3];
    
    int nonCommonVertex(int it1, int it2);

private:
    void read_mesh(char *fname);

    void construct_normal(TetraType *t, TriType tri_loc, int ittr_loc);

    int addVertex(istringstream &l,int halo=-1);
    int addVertex(double p[3],     int halo=-1);
    int addTetra(istringstream &l, int halo=-1);
    int addTetra(int p[4],         int halo=-1);
    int addEdges(int it, TetraType *t);
    int addTris(int it, TetraType *t);

    void fixNeighbors();
    
    int isTetra(int p[4], int start=0, int end=-1);
    int isTriangle(int v[3]);
    int isEdge(int v0, int v1);
    int isVertex(double p[3], int start=0, int end=-1);
    
    void ttrPeriodic(int it, int iv, int idim, int *ttr_new);
    void ttrPeriodic(int it, int iv2, bool per[3], int *ttr_new);

    void hashBoxesInit(int nVrt_target=1);
    void hashBoxesSort();
    void hashBoxAddVrt(int iv);
    void hashBoxAddTtr(int it);
    void hashBoxIndices(double *p, double *dx, int *i); 
    HashBox * hashBoxFromPoint(double *p, double *dx);

    double *** QRSplines_a[MAX_SPLINES];
    double **  QRSplines_b[MAX_SPLINES];
    double **  QRSplines_c[MAX_SPLINES];
    double **  QRSplines_d[MAX_SPLINES];
    double *   QRSplines_bare[MAX_SPLINES];
    
    int TtrSplinesND(int sType);
    
    bool haloDefined; 
    bool periodic[3];
    double *domain;
    
    double *hulXYZ; 
    
    int wrkSize; 
    void *wrk;
    vector<class Path> pth; 
    
    friend class Path;
};


class Path {

public:
    Path(); 
    Path(int *p, int n, double tau, Triangulation *t);

    int len();

    
private:
    vector<int> vrt;
    vector<int> edg; 
    friend class Triangulation;
    double tau; 
    
}; 
#endif // ifndef triangulate_h
