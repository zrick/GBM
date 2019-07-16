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
    
    std::cout << "Enter name of mesh file containing triangulation:";
    std::cin >> fname;
    
    read_mesh(fname);
}

Triangulation::Triangulation(char *fname){
    read_mesh(fname);
}

void Triangulation::read_mesh(char *fname){
    
    using namespace std; 
    // LOCAL DECLARATIONS
    int iLine=0,nEdgCommon,newEdges,newTris;
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
        lss >> p_vrt->p[0] >> p_vrt->p[1] >> p_vrt->p[2] >> p_vrt->e;
    }
    
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
        newEdges=add_edges(iTtr,p_ttr);  // update list of edges
        newTris =add_tris(iTtr,p_ttr);   // update list of triangles
        
    }
    
    // Fix Neighbors
    for ( int iTtr=0; iTtr<nTtr; ++iTtr ) {
        p_ttr=&ttr[iTtr];
        /*
        cout << "Tetra " << iTetra << " EDGES: ";
        for (int iEdg=0; iEdg<6; ++iEdg) cout << t->e[iEdg] << " ";
        cout << '\n';
         */
        
        for (int iAdj=0; iAdj <p_ttr->n_adjacent; iAdj++ )
        {
            int i_other = p_ttr->a[iAdj];
            /*
            cout << "ADJACENT TETRA:" << i_other <<  "Edges: ";
            for ( int iEdg=0; iEdg<6; ++iEdg) cout << tetras[i_other].e[iEdg] << " ";
            */
            nEdgCommon=0;
            for ( int iEdg=0; iEdg<6; ++ iEdg){
                for ( int iEdg2=0; iEdg2<6; ++iEdg2) {
                    if ( p_ttr->edg[iEdg] == ttr[i_other].edg[iEdg2] )
                        nEdgCommon++;
                }
            }
            
            if ( nEdgCommon == 3) {
                p_ttr->n[p_ttr->n_neighbor]=i_other;
                p_ttr->n_neighbor++;
                
                if ( p_ttr->n_neighbor > 4 ) {
                    cout << "ERROR: Internal Error; the number of neighbors cannot be larger than 4 \n \
                    Found " << p_ttr->n_neighbor << '\n';
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
    
    // Construct Hull; fix boudary property of vertices, edges, and terahedras
    for ( int iTri=0; iTri<nTri; ++iTri)
        if ( tri[iTri].bdy == true )
        {
            hul.push_back( &tri[iTri] );
          //  todo : fix according properties of vertices, edges and tetrahedras
            nHul++;
        }
    return;
}

int Triangulation::add_tris(int it, TetraType *t){
    int inew=0, v[3];
    int com[4][3] = { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}}; // combinations of vertices
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
                tri_loc.vrt[d]=v[d];                // add vertices
                tri_loc.c[d]+= vrt[v[d]].p[d] /3.;  // calculate centre
                if ( vrt[v[d]].nVrtTri >= MAX_VRT_TRI ) {
                    std::cout << "ERROR: Too many triangles (" <<vrt[v[d]].nVrtTri <<") for Vertex" << v[d] << '\n';
                    exit(EXIT_FAILURE);
                }
                if ( vrt[v[d]].nVrtTri > MAX_VRT_TRI ) {
                    std::cout << "ERROR: Too many triangles (" << vrt[v[d]].nVrtTri << ") on Vertex " << v[d] << '\n';
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
                std::cout << "ERROR: Too many Triangles (" << edg[tri[itri_loc].edg[d]].nEdgTri << ") on Edge " << tri[itri_loc].edg[d];
                exit(EXIT_FAILURE);
            }
            edg[tri[itri_loc].edg[d]].tri[edg[tri[itri_loc].edg[d]].nEdgTri]=itri_loc;
            edg[tri[itri_loc].edg[d]].nEdgTri++;
        }
    }
    return inew;
}

int Triangulation::is_triangle(int v[3]){
    for (int i=0; i<nTri; ++i)
        if ( tri[i].vrt[0] == v[0] && tri[i].vrt[1] == v[1] && tri[i].vrt[2] == v[2] )
            return(i);
    return (-1);
}

int Triangulation::add_edges(int it, TetraType *t){
    
    using namespace std; 
    int v0,v1, inew=0,iEdg,iEdg_loc=0,iTtr_loc,i2;
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
        
                    for (i2=0; i2<t->n_adjacent; ++i2){
                        if ( t->a[i2] == iTtr_loc ) break;
                    }
                    if (i2 == t->n_adjacent) {
                        if ( t->n_adjacent >= MAX_TTR_ADJ ) {
                            std::cout << "ERROR: Too many adjacent tetras (" << t->n_adjacent << ") on Tetra " << iTtr_loc;
                            exit(EXIT_FAILURE);
                        }
                        t->a[t->n_adjacent]=iTtr_loc;
                        t->n_adjacent++;
                    }
                   
                    t_other = &ttr[iTtr_loc];
           
                    for (i2=0;i2<t_other->n_adjacent; ++i2){
                        if ( t_other->a[i2] == it ) break;
                    }
                    if ( i2 == t_other->n_adjacent ) {
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
                for ( i2 =0; i2<v->nVrtEdg; ++i2) {
                    if ( v->edg[i2] == nEdg ) break;
                }
                if ( i2 == v->nVrtEdg) {
                    if ( v->nVrtEdg >= MAX_VRT_EDG ) {
                        cout << "ERROR: Too many Edges (" << v->nVrtEdg << ") on Vertex " << v0 << '\n';
                        exit(EXIT_FAILURE);
                    }
                    v->edg[v->nVrtEdg] = nEdg;
                    v->nVrtEdg++;
                }
               
                v=&vrt[v1];
                for ( i2=0; i2<v->nVrtEdg; ++i2){
                    if ( v->edg[i2] == nEdg) break;
                }
                if ( i2 == v->nVrtEdg ){
                    if ( v->nVrtEdg >= MAX_VRT_EDG ) {
                        cout << "ERROR: Too many Edges (" << v->nVrtEdg << ") on Vertex " << v1 << '\n';
                        exit(EXIT_FAILURE);
                    }
                    v->edg[v->nVrtEdg] = nEdg;
                    v->nVrtEdg++;
                }
                /*if (find (v->edges.begin(),v->edges.end(),nEdg-1) == v->edges.end() ) {
                    v->edges.push_back(nEdg-1);
                    v->nEdg++;
                }*/
                inew++;
                nEdg++;
            }
            iEdg_loc++;
        }
    }
    return inew;
}

int Triangulation::is_edge(int v0, int v1){
    
    for ( int iEdg=0; iEdg < nEdg; ++iEdg ) {
        if ( edg[iEdg].vrt[0] == v0 && edg[iEdg].vrt[1] == v1 )
            return(iEdg);
    }
    // 
    return ( -1 );
}

void Triangulation::printGrid(){
    
    using namespace std;
    VertexType *v;
    TetraType *t;

    cout << '\n' << nVrt << " VERTICES\n==============\n";
    for ( int i=0; i<nVrt; ++i){
        v=&vrt[i];
        cout << "Vrt " << i << "(";
        for (int i2=0; i2<3; ++i2) cout << v->p[i2] <<";";
        cout <<  "):" << v->nVrtEdg << "Edg: " ;
        for (int i2=0; i2<v->nVrtEdg; ++i2) cout << v->edg[i2] <<";";
        cout << v->nVrtTri << "Tri: ";
        for (int i2=0; i2<v->nVrtTri; ++i2) cout << v->tri[i2] <<";";
        cout << v->nVrtTtr << "Ttr: ";
        for (int i2=0; i2<v->nVrtTtr; ++i2) cout << v->ttr[i2] <<";";
        cout << '\n' ;
    }
    
    cout << '\n' << nEdg << " EDGES\n==============\n";
    for ( int i=0; i<nEdg; ++i){
        cout << "Edg " << i <<"("<< edg[i].vrt[0] << '-' <<  edg[i].vrt[1] << "):"<<  edg[i].nEdgTri << "Tri: ";
        for    (int n=0; n<edg[i].nEdgTri;   ++n) {
            cout << edg[i].tri[n] << " ";
        }
        cout << ";  " << edg[i].nEdgTtr << "Ttr: " ;
        for    (int n=0; n<edg[i].nEdgTtr; ++n) {
            cout << edg[i].ttr[n] << " ";
        }
        cout << '\n';
    }
    
    cout << '\n' << nTri << " TRIANGLES\n=================\n";
    for ( int i=0; i<nTri; ++i){
        cout <<"Tri " << i <<":(";
        for ( int d=0; d<3; ++d ) cout << tri[i].vrt[d] << ( d<2? "-" : ""  );
        cout << ") Edg: ";
        for ( int d=0; d<3; ++d ) cout << tri[i].edg[d] << ";";
        cout << " Ttr: ";
        for (int d=0; d<2; ++d) cout << tri[i].ttr[d] << ";";
        cout <<'\n';
    }

    
    cout << '\n' << nTtr <<" TETRAHEDRONS\n=================\n";
    for ( int i=0; i<nTtr; ++i){
        t=&ttr[i];
        cout << "Ttr " << i <<"(";
        for (int i2=0; i2<3; ++i2) cout << t->c[i2] << " ";
        cout <<") Vrt:";
        for (int i2=0; i2<4; ++i2) cout << t->vrt[i2] << ";";
        cout << "Edg: ";
        for (int i2=0; i2<6; ++i2) cout << t->edg[i2] << ";";
        cout << "Tri: ";
        for (int i2=0; i2<4; ++i2) cout << t->tri[i2] << ";";
        cout << "Adj: ";
        for ( int i2=0;i2<t->n_adjacent; ++i2) cout << t->a[i2] << ";";
        cout << "Ngh: ";
        for ( int i2=0;i2<t->n_neighbor; ++i2) cout << t->n[i2] << ";";
        cout << '\n';
    }
    
    
}


void Triangulation::sort4(int a[4])
{
    using namespace std;
    
    if (a[0]>a[1]) swap(a[0],a[1]);  // first two values in order
    if (a[1]>a[2]) {
        swap(a[1],a[2]);
        if (a[0]>a[1]) swap(a[0],a[1]);
    }                                // first three values in order
    if (a[2]>a[3]) {
        swap(a[2],a[3]);
        if (a[1]>a[2]) swap(a[1],a[2]);
        if (a[0]>a[1]) swap(a[0],a[1]);
    }                                // all (four) values in order
    
    return;
}
