//
//  types.h
//  GBM
//
//  Created by Cedrick Ansorge on 30.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef types_h
#define types_h

#include <stdio.h>
#include <cstring>

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
    bool hasAttribute(string &n);

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
    Namelist(string &fname);
    
    bool   hasVal(string group, string var);
    
    int    getVal_int(string group,string name);
    double getVal_dbl(string group,string name);
    string getVal_str(string group,string name);
    bool   getVal_bool(string group,string name);
    
    void getList_int(string group, string var,int    *list, const int nmax);
    void getList_dbl(string group, string var,double *list, const int nmax);

    
private:
    //Functions
    void read_namelist(string s);
    string getValStr(string group,string var);
    
    // Members
    std::vector<Group> grp;
    std::vector<string>grpNames;

};

typedef struct VertexType{
    int idx;
    double p[3];   // position
    double e;      // elevation on unit-paraboloid ( 2-norm of p; )
    int atl;       // Atlas label
    double vol,dia;// Volume and diameter
    int edg[MAX_VRT_EDG];
    int tri[MAX_VRT_TRI];
    int ttr[MAX_VRT_TTR];
    int nVrtEdg, nVrtTri, nVrtTtr;
    bool bdy;
} VertexType;

typedef struct EdgeType{
    int idx;
    int vrt[2];      // Vertices connected
    int ttr[MAX_EDG_TTR];      // Tetrahedrons sharing this edge
    int tri[MAX_EDG_TRI];     // Triangles sharing this edge
    int nEdgTtr, nEdgTri;
    bool bdy;
} EdgeType;

typedef struct TriType{
    int idx;
    int vrt[3];    // Vertices
    int edg[3];    // Edges
    int ttr[2];    // Tetrahedrons
    double c[3];   // Centroid
    double n[2][3];// Outer normals w.r.t. Tetrahedra ttr[0] and ttr[1]
    double area; 
    bool bdy;      // true if part of external surface
} TriType;

typedef struct TetraType{
    int idx;
    int atl;
    int vrt[4]; // Vertices
    int tri[4]; // Triangles
    int edg[6]; // Edges
    double c[3];// Centroid point
    double vol; // Volume
    int n[4];   // Neighbors (sharing a surface)
    int a[MAX_TTR_ADJ];  // adjacents (sharing an edge)
    int n_neighbor;
    int n_adjacent;
    bool bdy;
} TetraType;



#endif /* types_h */
