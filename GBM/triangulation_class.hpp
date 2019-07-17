//
//  triangulation_class.hpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef triangulation_class_hpp

#include "main.hpp"
#include "utils.hpp"

using namespace std; 

#define MAX_VRT_TTR 128
#define MAX_VRT_TRI 128
#define MAX_VRT_EDG 128

#define MAX_EDG_TRI 32
#define MAX_EDG_TTR 32

#define MAX_TTR_ADJ 64

typedef struct VertexType{
    double p[3];   // position
    double e;      // elevation on unit-paraboloid ( 2-norm of p; )
    int edg[MAX_VRT_EDG];
    int tri[MAX_VRT_TRI];
    int ttr[MAX_VRT_TTR];
    int nVrtEdg, nVrtTri, nVrtTtr;
    bool bdy;
} VertexType;

typedef struct EdgeType{
    int vrt[2];      // Vertices connected
    int ttr[MAX_EDG_TTR];      // Tetrahedrons sharing this edge
    int tri[MAX_EDG_TRI];     // Triangles sharing this edge
    int nEdgTtr, nEdgTri;
    bool bdy;
} EdgeType;

typedef struct TriType{
    int vrt[3];    // Vertices
    int edg[3];    // Edges
    int ttr[2];    // Tetrahedrons
    double c[3];   // Centroid
    double n[3];   // Normal
    bool bdy;      // true if part of external surface
} TriType;

typedef struct TetraType{
    int vrt[4];     // Vertices
    int tri[4];     // Triangles
    int edg[6];     // Edges
    double c[3];  // Centroid point
    int n[4];   // Neighbors (sharing a surface)
    int a[MAX_TTR_ADJ];  // adjacents (sharing an edge)
    int n_neighbor;
    int n_adjacent;
    bool bdy;
} TetraType;

class Triangulation{
    
public:
    // Constructors
    Triangulation();
    Triangulation(char *fname);
 
    void printGrid(int level);
    
    void printVertex(int i, int level);
    void printEdge(int i, int level);
    void printTriangle(int i, int level);
    void printTetrahedron(int i, int level);
    
    TetraType *ttr;      // 3D Tetrahedrons
    VertexType *vrt;     // 0D Vertex Points
    std::vector<TriType>    tri;     // 2D Triangles (surfaces of Tetrahedrons)
    std::vector<EdgeType>   edg;    // 1D Edges of Tetrahedrons and Triangles

    std::vector<TriType *>  hul;     // Exterior Surface of Triangulation
    // Members
    int nDim, nVrt, nEdg, nTri, nTtr, nHul;
    
private:
    void read_mesh(char *fname);
    int add_edges(int it, TetraType *t);
    int add_tris(int it, TetraType *t);
    int is_triangle(int v[3]);
    int is_edge(int v0, int v1);
    
};

#define triangulation_class_hpp

#endif /* triangulation_class_hpp */
