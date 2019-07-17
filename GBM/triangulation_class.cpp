//
//  triangulation_class.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//
#include "triangulation_class.hpp"

Triangulation::Triangulation(){
    char fname[MAX_CHAR_LEN];
    
    cout << "Enter name of mesh file containing triangulation:";
    cin >> fname;
    
    read_mesh(fname);
}

Triangulation::Triangulation(char *fname){
    read_mesh(fname);
}

void Triangulation::read_mesh(char *fname){
    // LOCAL DECLARATIONS
    int iLine=0,nEdgCommon;
    struct stat buffer;
    string line;
    ifstream tfile(fname);
    VertexType *p_vrt;
    TetraType *p_ttr;
    
    // check if file exists
    if ( stat (fname, &buffer) != 0) {
        cout<<"Error: file \'" << &fname[0] << "\' not found \n";
    } else {
        cout << "Reading mesh from file \'" << &fname[0] << "\'\n";
    }

    
    // READ VERTICES
    {
        getline(tfile,line); iLine++;
        istringstream lss(line);
        lss >> nDim; nDim--;   // decrement by one, because it includes elevation
        if ( nDim != 3 ) {
            cout << "ERROR: in " << fname << ":" << iLine+1 << '\n';
            cout << "       Dimensions of " << nDim << " not supported" << '\n';
            exit(EXIT_FAILURE);
        }
    }
    
    {
        getline(tfile,line); iLine++;
        istringstream lss(line);
        lss >> nVrt;
    }
  
    cout << "Dimensions: " << nDim << '\n';
    cout << "Reading " << nVrt <<"Vertices \n";
    vrt = (VertexType *) malloc(nVrt*sizeof(VertexType));
    for ( int iVrt=0; iVrt<nVrt; ++iVrt) {
        p_vrt=&vrt[iVrt];
        getline(tfile,line);
        istringstream lss(line);
        p_vrt->nVrtEdg=0;
        p_vrt->nVrtTri=0;
        p_vrt->nVrtTtr=0;
        lss >> p_vrt->p[0] >> p_vrt->p[1] >> p_vrt->p[2] >> p_vrt->e;\
        p_vrt->nVrtEdg=0;
        p_vrt->nVrtTri=0;
        p_vrt->nVrtTtr=0;
    }
    
    cout << "Processed Vertices\n";
    // READ TETRAHEDRONS
    {
        getline(tfile,line); iLine++;
        istringstream lss(line);
        lss >> nTtr;
    }
    
    cout << "Reading " << nTtr << "Tetrahedrons" << '\n';
    ttr = (TetraType *) malloc (nTtr * sizeof(TetraType));
    
    
    for ( int iTtr=0; iTtr<nTtr; ++iTtr) {
        p_ttr=&ttr[iTtr];
        getline(tfile,line);
        istringstream lss(line);
        lss >> p_ttr->vrt[0] >> p_ttr->vrt[1] >> p_ttr->vrt[2] >> p_ttr->vrt[3];
        sort4(p_ttr->vrt);
        p_ttr->n_adjacent=0;
        p_ttr->n_neighbor=0;
        
        // Calculate Centroid of tetrahedron
        for ( int idim=0; idim<3; ++idim) {
            p_ttr->c[idim]=0.;
            for ( int ivert=0; ivert<4; ++ivert) {
                p_ttr->c[idim] += vrt[p_ttr->vrt[ivert]].p[idim];
            }
            p_ttr->c[idim] *= 0.25;
        }
        
        for ( int i=0; i<4; ++i) {
            p_vrt=&vrt[p_ttr->vrt[i]];
            if ( p_vrt->nVrtTtr > MAX_VRT_TTR )  {
                cout << "ERROR: Too many tetrahedrons (" <<p_vrt->nVrtTtr <<") on Vertex" << p_ttr->vrt[i] << '\n';
                exit(EXIT_FAILURE);
            }
            p_vrt->ttr[p_vrt->nVrtTtr]=iTtr;
            p_vrt->nVrtTtr++;
        }
        add_edges(iTtr,p_ttr);  // update list of edges
        add_tris(iTtr,p_ttr);   // update list of triangles
        
    }
    cout << "Processed Tetrahedrons \n";
    // Fix Neighbors
    for ( int iTtr=0; iTtr<nTtr; ++iTtr ) {
        p_ttr=&ttr[iTtr];

        for (int iAdj=0; iAdj < p_ttr->n_adjacent; iAdj++ )
        {
            int i_other = p_ttr->a[iAdj];
            nEdgCommon=0;
            for ( int iEdg=0; iEdg<6; ++ iEdg)   // count all common edges
                if ( ifind(ttr[i_other].edg,6,p_ttr->edg[iEdg]) != 6 )
                    nEdgCommon++;
            
            if ( nEdgCommon == 3) {
                p_ttr->n[p_ttr->n_neighbor]=i_other;
                p_ttr->n_neighbor++;
                
                if ( p_ttr->n_neighbor > 4 ) {
                    cout << "ERROR: Internal Error; the number of neighbors cannot be larger than 4 \n \
                    Found " << p_ttr->n_neighbor << "for Tetra " << iTtr << '\n';
                    exit(EXIT_FAILURE);
                }
            } else if ( nEdgCommon != 1 ) {
                cout << "ERROR: Internal Error the number of common Edges among Tetrahedrons \n \
                Must either be 1 or 3;\n       Found " << nEdgCommon;
                 exit(EXIT_FAILURE);
            }
        }
        if ( p_ttr-> n_neighbor <1 || p_ttr->n_neighbor > 4){
            cout << "ERROR: Internal Error the number of neighbors for a Tetrahedron must be 2,3, or 4; found " << p_ttr->n_neighbor << "\n";
        }
    }
    
    cout << "Fixed Neighbors\n";
    
    // Construct Hull; fix boudary property of vertices, edges, and terahedras
    for ( int iTri=0; iTri<nTri; ++iTri)
        if ( tri[iTri].bdy == true )
        {
            hul.push_back( &tri[iTri] );
          //  TODO : fix according properties of vertices, edges and tetrahedras
            nHul++;
        }

    cout << "Calculated Hull\n";
    return;
}

int Triangulation::add_tris(int it, TetraType *t){
    int inew=0,v[3];
    int com[4][3] = { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}}; // 4 possible combinations of the three vertices
    int itri_loc;
  
    for ( int i=0; i<4; ++i){
        for (int d=0; d<3; ++d) v[d]=t->vrt[com[i][d]];
        itri_loc=is_triangle(v);
        
        if ( itri_loc < 0) { // Triangle needs to be added
            TriType tri_loc;
            
            itri_loc=nTri;
            tri_loc.edg[0]=is_edge(v[0],v[1]);
            tri_loc.edg[1]=is_edge(v[0],v[2]);
            tri_loc.edg[2]=is_edge(v[1],v[2]);

            for (int d=0; d<3; ++d) {
                tri_loc.vrt[d]=v[d];         // add vertices
                tri_loc.c[d]=0.;                     // calculate centre
                for(int dd=0;dd<3;++dd)
                    tri_loc.c[d] += vrt[v[dd]].p[d] /3.;
                if ( vrt[v[d]].nVrtTri >= MAX_VRT_TRI ) {
                    cout << "ERROR: Too many triangles (" <<vrt[v[d]].nVrtTri <<") for Vertex" << v[d] << '\n';
                    exit(EXIT_FAILURE);
                }
                if ( vrt[v[d]].nVrtTri > MAX_VRT_TRI ) {
                    cout << "ERROR: Too many triangles (" << vrt[v[d]].nVrtTri << ") on Vertex " << v[d] << '\n';
                    exit(EXIT_FAILURE);
                }
                vrt[v[d]].tri[vrt[v[d]].nVrtTri]=itri_loc;
                vrt[v[d]].nVrtTri++;
            }
            
            tri_loc.ttr[0]=it;
            tri_loc.ttr[1]=-1;
            tri_loc.bdy=true;
            tri.push_back(tri_loc);
            
            nTri++;
            inew++;
        }else{
            tri[itri_loc].ttr[1]=it;
            tri[itri_loc].bdy=false;
        }
        
        t->tri[i]=itri_loc;
        for (int d=0; d<3; ++d) {
            if( edg[tri[itri_loc].edg[d]].nEdgTri >= MAX_EDG_TRI ) {
                cout << "ERROR: Too many Triangles (" << edg[tri[itri_loc].edg[d]].nEdgTri << ") on Edge " << tri[itri_loc].edg[d];
                exit(EXIT_FAILURE);
            }
            edg[tri[itri_loc].edg[d]].tri[edg[tri[itri_loc].edg[d]].nEdgTri]=itri_loc;
            edg[tri[itri_loc].edg[d]].nEdgTri++;
        }
    }
    return inew;
}

int Triangulation::is_triangle(int v[3]){
    TriType *t_loc;
    VertexType *v_loc;
    for ( int i=0; i<3; ++i){
        v_loc=&vrt[v[i]];
        for ( int i2=0; i2 < v_loc->nVrtTri; ++i2)
        {
            t_loc=&tri[v_loc->tri[i2]];
            if ( t_loc->vrt[0] == v[0] && t_loc->vrt[1] == v[1] && t_loc->vrt[2] == v[2] )
                return v_loc->tri[i2];
        }
    }
    return -1;
    
  /*
    for (int i=0; i<nTri; ++i)
        if ( tri[i].vrt[0] == v[0] && tri[i].vrt[1] == v[1] && tri[i].vrt[2] == v[2] )
            return(i);
    return (-1);
   */
}

int Triangulation::add_edges(int it, TetraType *t){
    int v0,v1, inew=0,iEdg,iEdg_loc=0,iTtr_loc;
    TetraType *t_other;
    VertexType *v;
    
    iEdg_loc=0;
    for ( int i0=0; i0<4; ++i0){
        for ( int i1=i0+1; i1<4; ++i1) {
            v0=t->vrt[i0];
            v1=t->vrt[i1];
            if ( v0 > v1 ) swap(v0,v1);
       
            iEdg = is_edge(v0,v1);
            if ( iEdg >= 0 )      // Edge exists already; add information about Tetrahedron
            {
                // recover other tetrahedrons from edge properties and tell them we are adjacent
                for ( int it_loc=0; it_loc<edg[iEdg].nEdgTtr; ++it_loc) {
                    iTtr_loc=edg[iEdg].ttr[it_loc];
                    if (ifind(t->a,t->n_adjacent,iTtr_loc) == t->n_adjacent) {
                        if ( t->n_adjacent >= MAX_TTR_ADJ ) {
                            cout << "ERROR: Too many adjacent tetras (" << t->n_adjacent << ") on Tetra " << iTtr_loc;
                            exit(EXIT_FAILURE);
                        }
                        t->a[t->n_adjacent]=iTtr_loc;
                        t->n_adjacent++;
                    }
                   
                    t_other = &ttr[iTtr_loc];
           
                    if ( ifind(t_other->a,t_other->n_adjacent,it) == t_other->n_adjacent ) {
                        if ( t_other->n_adjacent >= MAX_TTR_ADJ ) {
                            cout << "ERROR: Too many adjacent tetras (" << t_other->n_adjacent << ") on Tetra " << iTtr_loc;
                            exit(EXIT_FAILURE);
                        }
                        t_other->a[t_other->n_adjacent]=it;
                        t_other->n_adjacent++;
                    }
                }
                if (edg[iEdg].nEdgTtr >= MAX_EDG_TTR ) {
                    cout << "ERROR: Too many Tetras (" << edg[iEdg].nEdgTtr << ") on Edge " << iEdg << '\n';
                    exit(EXIT_FAILURE);
                }
                edg[iEdg].ttr[edg[iEdg].nEdgTtr]=it;
                edg[iEdg].nEdgTtr++;
                t->edg[iEdg_loc]=iEdg;
            }
            else                 // Edge does not exist yet; needs to be created
            {
                EdgeType edg_loc;
                
                edg_loc.vrt[0]=v0;
                edg_loc.vrt[1]=v1;
                edg_loc.ttr[0]=it;
                edg_loc.nEdgTtr=1;
                edg_loc.nEdgTri=0;
                edg.push_back(edg_loc);
                t->edg[iEdg_loc]=nEdg;
                
                v=&vrt[v0];
                if ( ifind(v->edg,v->nVrtEdg,nEdg) == v->nVrtEdg) {
                    if ( v->nVrtEdg >= MAX_VRT_EDG ) {
                        cout << "ERROR: Too many Edges (" << v->nVrtEdg << ") on Vertex " << v0 << '\n';
                        exit(EXIT_FAILURE);
                    }
                    v->edg[v->nVrtEdg] = nEdg;
                    v->nVrtEdg++;
                }
               
                v=&vrt[v1];
                if ( ifind(v->edg,v->nVrtEdg,nEdg) == v->nVrtEdg ){
                    if ( v->nVrtEdg >= MAX_VRT_EDG ) {
                        cout << "ERROR: Too many Edges (" << v->nVrtEdg << ") on Vertex " << v1 << '\n';
                        exit(EXIT_FAILURE);
                    }
                    v->edg[v->nVrtEdg] = nEdg;
                    v->nVrtEdg++;
                }
                inew++;
                nEdg++;
            }
            iEdg_loc++;
        }
    }
    return inew;
}

int Triangulation::is_edge(int v0, int v1){
    // check if one of the vertices already holds this edge
    for ( int iEdg=0; iEdg < vrt[v0].nVrtEdg; ++iEdg )
        if ( edg[vrt[v0].edg[iEdg]].vrt[1] == v1 )
            return vrt[v0].edg[iEdg];
    return ( -1 );
}

void Triangulation::printGrid(int level ){
    /* PARAMETER level:
        0 - bulk grid information;
        1 - basic connectivity;
        2 - detailed connectivity   */

    if ( level >=0 ) cout << '\n' << nVrt << " VERTICES\n==============\n";
    if ( level > 0 ) for ( int i=0; i<nVrt; ++i) printVertex(i,level);
    
    if ( level >=0 )cout << nEdg << " EDGES\n==============\n";
    if ( level > 0 ) for ( int i=0; i<nEdg; ++i) printEdge(i,level);
    
    if ( level >=0 ) cout << nTri << " TRIANGLES\n=================\n";
    if ( level > 0 ) for ( int i=0; i<nTri; ++i) printTriangle(i,level);
    
    if ( level >=0 ) cout << nTtr <<" TETRAHEDRONS\n=================\n";
    if ( level > 0 ) for ( int i=0; i<nTtr; ++i) printTetrahedron(i,level);

}

void Triangulation::printVertex(int iv, int level) {
    VertexType *v=&vrt[iv];
    
    cout << "Vrt " << iv << "(";
    if (level >0 ) for (int i2=0; i2<3; ++i2) cout << v->p[i2] <<";";
    cout <<  "):" << v->nVrtEdg << "Edg: " ;
    if ( level>1) for (int i2=0; i2<v->nVrtEdg; ++i2) cout << v->edg[i2] <<";";
    cout << v->nVrtTri << "Tri: ";
    if (level>1)  for (int i2=0; i2<v->nVrtTri; ++i2) cout << v->tri[i2] <<";";
    cout << v->nVrtTtr << "Ttr: ";
    if ( level>1) for (int i2=0; i2<v->nVrtTtr; ++i2) cout << v->ttr[i2] <<";";
    cout << "\n" ;
}
void Triangulation::printEdge(int i, int level) {
    EdgeType *e=&edg[i];
    
    cout << "Edg " << i <<"("<< e->vrt[0] << '-' <<  e->vrt[1] << "):"<<  e->nEdgTri << "Tri: ";
    if ( level > 1 ) for    (int n=0; n<e->nEdgTri;   ++n) cout << e->tri[n] << " ";
    cout << ";  " << edg[i].nEdgTtr << "Ttr: " ;
    if ( level > 1 ) for    (int n=0; n<e->nEdgTtr; ++n) cout << e->ttr[n] << " ";
    cout << "\n";
}

void Triangulation::printTriangle(int i, int level) {
    TriType*t=&tri[i];
    
    cout <<"Tri " << i <<":(";
    if ( level >0 ) for ( int d=0; d<3; ++d ) cout << t->vrt[d] << ( d<2? "-" : ")"  );
    if ( level >1 ) {
        cout << " Edg: ";
        for ( int d=0; d<3; ++d ) cout << t->edg[d] << ";";
        cout << " Ttr: ";
        for (int d=0; d<2; ++d) cout << t->ttr[d] << ";";
    }
    cout << '\n';
}

void Triangulation::printTetrahedron(int i, int level){
    TetraType *t=&ttr[i];

    cout << "Ttr " << i <<"(";
    for (int i2=0; i2<3; ++i2) cout << t->c[i2] << " ";
    cout <<") Vrt:";
    for (int i2=0; i2<4; ++i2) cout << t->vrt[i2] << ";";
    if ( level > 1 ) {
        cout << "Edg: ";
        if ( level> 1 ) for (int i2=0; i2<6; ++i2) cout << t->edg[i2] << ";";
        cout << "Tri: ";
        if ( level> 1 ) for (int i2=0; i2<4; ++i2) cout << t->tri[i2] << ";";
    }
    cout << "(" << t->n_adjacent << " Adj): ";
    if ( level> 1 ) for ( int i2=0;i2<t->n_adjacent; ++i2) cout << t->a[i2] << ";";
    cout << "(" << t->n_neighbor << " Ngh): ";
    if ( level> 1 ) for ( int i2=0;i2<t->n_neighbor; ++i2) cout << t->n[i2] << ";";
    cout << "\n";
}

