//
//  triangulate.cpp
//  triangulate
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

//
//  triangulate.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//
#include "triangulate.hpp"

Triangulation::Triangulation(){
    /*char fname[MAX_CHAR_LEN];
    
    cout << "Enter name of mesh file containing triangulation:";
    cin >> fname;
    
    read_mesh(fname);
     */
    return; 
}

Triangulation::Triangulation(char *fname){
    read_mesh(fname);
}

void Triangulation::read_mesh(char *fname){
    // LOCAL DECLARATIONS
    int iLine=0,nEdgCommon;
    string line;
    VertexType *p_vrt;
    TetraType *p_ttr;
    TriType *p_tri;
    
    // check if file exists
    if ( file_exist(fname) )
        cout << "Reading mesh from file \'" << &fname[0] << "\'\n";
    
    ifstream tfile(fname);
    
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
    for ( int i=0; i<3; ++i ){
        box[2*i]= 9e9;
        box[2*i+1]=-9e9;
    }
    
    for ( int iVrt=0; iVrt<nVrt; ++iVrt) {
        p_vrt=&vrt[iVrt];
        p_vrt->idx=iVrt;
        getline(tfile,line);
        istringstream lss(line);
        p_vrt->nVrtEdg=0;
        p_vrt->nVrtTri=0;
        p_vrt->nVrtTtr=0;
        lss >> p_vrt->p[0] >> p_vrt->p[1] >> p_vrt->p[2] >> p_vrt->e;\
        p_vrt->nVrtEdg=0;
        p_vrt->nVrtTri=0;
        p_vrt->nVrtTtr=0;
        
        box[0] = p_vrt->p[0] < box[0]? p_vrt->p[0] : box[0];
        box[1] = p_vrt->p[0] > box[1]? p_vrt->p[0] : box[1];
        box[2] = p_vrt->p[1] < box[2]? p_vrt->p[1] : box[2];
        box[3] = p_vrt->p[1] > box[3]? p_vrt->p[1] : box[3];
        box[4] = p_vrt->p[2] < box[4]? p_vrt->p[2] : box[4];
        box[5] = p_vrt->p[2] > box[5]? p_vrt->p[2] : box[5];
    
    }
    
    // READ TETRAHEDRONS
    {
        getline(tfile,line); iLine++;
        istringstream lss(line);
        lss >> nTtr;
    }
    
    cout << "Reading " << nTtr << "Tetrahedrons" << '\n';
    ttr = (TetraType *) malloc (nTtr * sizeof(TetraType));
    
    volume = 0. ;
    for ( int iTtr=0; iTtr<nTtr; ++iTtr) {
        p_ttr=&ttr[iTtr];
        getline(tfile,line);
        istringstream lss(line);
        lss >> p_ttr->vrt[0] >> p_ttr->vrt[1] >> p_ttr->vrt[2] >> p_ttr->vrt[3];
        sort4(p_ttr->vrt);
        p_ttr->idx=iTtr;
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
        p_ttr->vol=tetraVolume(vrt[p_ttr->vrt[0]].p, vrt[p_ttr->vrt[1]].p,
                               vrt[p_ttr->vrt[2]].p, vrt[p_ttr->vrt[3]].p);
        volume += p_ttr->vol; 
        add_edges(iTtr,p_ttr);  // update list of edges
        add_tris(iTtr,p_ttr);   // update list of triangles
    }
    
    // Fix Neighbors
    cout << "Fixing Neighbors \n";
    for ( int iTtr=0; iTtr<nTtr; ++iTtr ) {
        p_ttr=&ttr[iTtr];
        
        for (int iAdj=0; iAdj < p_ttr->n_adjacent; iAdj++ )
        {
            int i_other = p_ttr->a[iAdj];
            nEdgCommon=0;
            for ( int iEdg=0; iEdg<6; ++ iEdg)   // count all common edges
                if ( afind(ttr[i_other].edg,6,p_ttr->edg[iEdg]) != 6 )
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
 
    // Construct Hull; fix boudary property of vertices, edges, and terahedras
    for ( int iTri=0; iTri<nTri; ++iTri)
        if ( tri[iTri].bdy == true )
        {
            p_tri=&tri[iTri];
            hul.push_back( p_tri );
            for (int d=0; d<3; ++d) {
                edg[p_tri->edg[d]].bdy=true;
                vrt[p_tri->vrt[d]].bdy=true;
            }
            ttr[p_tri->ttr[0]].bdy=true;
            if ( p_tri->ttr[1] > 0 ) {
                cout << "ERROR: inconsistent boundary property on Triangle" << iTri << " (" << p_tri->vrt[0] <<"-" <<p_tri->vrt[1] << "-" << p_tri->vrt[2] <<'\n';
                exit(EXIT_FAILURE);
            }
            nHul++;
        }
    
    cout << "Calculated Hull; contains " << nHul << " triangles. \n";
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
            tri_loc.idx=nTri;
            itri_loc=nTri;
            tri_loc.edg[0]=is_edge(v[0],v[1]);
            tri_loc.edg[1]=is_edge(v[0],v[2]);
            tri_loc.edg[2]=is_edge(v[1],v[2]);
            
            for (int d=0; d<3; ++d) {
                tri_loc.vrt[d]=v[d];                 // add vertices
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
            construct_normal(t,tri_loc,0);
            
            tri_loc.ttr[1]=-1;
            tri_loc.bdy=true;
            tri.push_back(tri_loc);
            
            nTri++;
            inew++;
        }else{
            tri[itri_loc].ttr[1]=it;
            construct_normal(t,tri[itri_loc],1);
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
                    if (afind(t->a,t->n_adjacent,iTtr_loc) == t->n_adjacent) {
                        if ( t->n_adjacent >= MAX_TTR_ADJ ) {
                            cout << "ERROR: Too many adjacent tetras (" << t->n_adjacent << ") on Tetra " << iTtr_loc;
                            exit(EXIT_FAILURE);
                        }
                        t->a[t->n_adjacent]=iTtr_loc;
                        t->n_adjacent++;
                    }
                    
                    t_other = &ttr[iTtr_loc];
                    
                    if ( afind(t_other->a,t_other->n_adjacent,it) == t_other->n_adjacent ) {
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
                edg_loc.idx=nEdg;
                edg_loc.vrt[0]=v0;
                edg_loc.vrt[1]=v1;
                edg_loc.ttr[0]=it;
                edg_loc.nEdgTtr=1;
                edg_loc.nEdgTri=0;
                edg.push_back(edg_loc);
                t->edg[iEdg_loc]=nEdg;
                
                v=&vrt[v0];
                if ( afind(v->edg,v->nVrtEdg,nEdg) == v->nVrtEdg) {
                    if ( v->nVrtEdg >= MAX_VRT_EDG ) {
                        cout << "ERROR: Too many Edges (" << v->nVrtEdg << ") on Vertex " << v0 << '\n';
                        exit(EXIT_FAILURE);
                    }
                    v->edg[v->nVrtEdg] = nEdg;
                    v->nVrtEdg++;
                }
                
                v=&vrt[v1];
                if ( afind(v->edg,v->nVrtEdg,nEdg) == v->nVrtEdg ){
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

void Triangulation::writeGrid(string grid_file, string grid_format) {

    /* PARAMETER format:
     0 - bulk grid information;
     1 - basic connectivity;
     2 - detailed connectivity
     3 - VTK XML file*/

    int format;
    double *data;
    int *idata;
    int iv,it,d,data_size;
    
    format=FMT_ASC_BASIC;
    if ( ! grid_format.compare("ASCII_BASIC"))
        format=FMT_ASC_BASIC;
    else if (! grid_format.compare("ASCII_FULL"))
        format=FMT_ASC_FULL;
    else if (! grid_format.compare("XML_VTK"))
        format=FMT_XML_VTK;

    cout << "WRITING GRID: " << grid_file << " "<< format <<" \n";
    
    if ( format <= FMT_ASC_FULL ) { // ASCII OUTPUT
        if ( format >=FMT_ASC_BASIC ) cout << '\n' << nVrt << " VERTICES\n==============\n";
        if ( format > FMT_ASC_BASIC ) for ( int i=0; i<nVrt; ++i) printVertex(&vrt[i],format);
    
        if ( format >=FMT_ASC_BASIC )cout << nEdg << " EDGES\n==============\n";
        if ( format > FMT_ASC_BASIC ) for ( int i=0; i<nEdg; ++i) printEdge(&edg[i],format);
    
        if ( format >=FMT_ASC_BASIC ) cout << nTri << " TRIANGLES\n=================\n";
        if ( format > FMT_ASC_BASIC ) for ( int i=0; i<nTri; ++i) printTriangle(&tri[i],format);
    
        if ( format >=FMT_ASC_BASIC ) cout << nTtr <<" TETRAHEDRONS\n=================\n";
        if ( format > FMT_ASC_BASIC ) for ( int i=0; i<nTtr; ++i) printTetrahedron(&ttr[i],format);
    }
    else if ( format == FMT_XML_VTK ) { // VTK XML FILE
        vector<string> att;
        ofstream gfile;
       
        data_size = 3*nVrt > 4*nTtr? 3*nVrt : 4*nTtr;
        
        
        data =(double*)malloc(data_size*sizeof(double));
        idata=(int*)   malloc(data_size*sizeof(int));
        
        vtkXMLFileOpen(grid_file,nVrt,nTtr,gfile);

        // GRID INFORMATION -- VERTICES
        gfile << "<Points>\n";
        // <DataArray type="Float32" NumberOfComponents="3" format="ascii">
        att.push_back("type");              att.push_back("Float32");
        att.push_back("NumberOfComponents");att.push_back("3");
        att.push_back("format");            att.push_back("ascii");
       
        for (iv=0;iv<nVrt;++iv)
            for (d=0;d<3;++d)
                data[3*iv+d]=vrt[iv].p[d];
        vtkXMLWriteDataArray(gfile, att, 3*nVrt, data);
        gfile << "</Points>\n";

        att.clear();
        
        gfile << "<Cells>\n";
        
        
        // <DataArray type="Int32" Name="connectivity" format="ascii">
        att.push_back("type");  att.push_back("Int32");
        att.push_back("Name");  att.push_back("connectivity");
        att.push_back("format");att.push_back("ascii");
        for (it=0; it<nTtr; ++it)
            for (d=0;d<4;++d)
                idata[4*it+d] = ttr[it].vrt[d];
        vtkXMLWriteDataArray(gfile, att, 4*nTtr, idata);
        att.clear();
        
        // <DataArray type="UInt8" Name="types" format="ascii">
        att.push_back("type");  att.push_back("UInt8");
        att.push_back("Name");  att.push_back("types");
        att.push_back("format");att.push_back("ascii");
        for (it=0; it<nTtr; ++it)
            idata[it] =10;
        vtkXMLWriteDataArray(gfile, att, nTtr, idata);
        att.clear();
        
        // <DataArray type="Int32" Name="offsets" format="ascii">
        att.push_back("type");  att.push_back("Int32");
        att.push_back("Name");  att.push_back("offsets");
        att.push_back("format");att.push_back("ascii");
        for (it=0; it<nTtr; ++it)
            idata[it] =(it+1)*4;
        vtkXMLWriteDataArray(gfile, att, nTtr, idata);
        att.clear();
        
        gfile << "</Cells>\n";
        
        
        gfile << "<PointData>\n";
 
        // <DataArray type="Int32" Name="Atlas_Region" format="ascii">
        att.push_back("type");  att.push_back("Int32");
        att.push_back("Name");  att.push_back("Atlas_Region");
        att.push_back("format");att.push_back("ascii");
        for(iv=0; iv<nVrt;++iv )
        {
            idata[iv]    =vrt[iv].atl;
            data[iv]     =vrt[iv].dia;
            data[nVrt+iv]=vrt[iv].vol;
            cout << iv << " " << vrt[iv].vol << '\n';
        }
        vtkXMLWriteDataArray(gfile,att,nVrt,idata);
        att[1]="Float32"; att[3]="Atlas_Diameter";
        vtkXMLWriteDataArray(gfile,att,nVrt,data);
        
        att[3]="Atlas_Volume";
        vtkXMLWriteDataArray(gfile, att, nVrt, &data[nVrt]);
        
        gfile << "</PointData>\n";
        
        gfile << "<CellData>\n";
        // <DataArray type="Int32" Name="Atlas_Region" format="ascii">
        att[1]="Int32"; att[3]="Atlas_Region";
        for ( it=0; it<nTtr; ++it) {
            idata[it]=ttr[it].atl;
            data[it] =ttr[it].vol;
        }
        vtkXMLWriteDataArray(gfile, att, nTtr, idata);
        
        att[1]="Float32";att[3]="Tetrahedron_Volume";
        vtkXMLWriteDataArray(gfile, att, nTtr, data);
            
        
        gfile << "</CellData>\n";

        // FOOTER of VTK XML FILE
        vtkXMLFileClose(gfile);
        
        att.clear();
        free(idata);
        free(data);
    }
}

void Triangulation::printHull(int level){
    if ( level > FMT_ASC_BASIC ) for ( int i=0; i<nVrt; ++i) if ( vrt[i].bdy ) printVertex(&vrt[i],level);
    
    if ( level >=FMT_ASC_BASIC )cout << " EDGES ON HULL\n==============\n";
    if ( level > FMT_ASC_BASIC ) for ( int i=0; i<nEdg; ++i) if ( edg[i].bdy ) printEdge(&edg[i],level);
    
    if ( level >=FMT_ASC_BASIC ) cout << " TRIANGLES ON HULL\n=================\n";
    if ( level > FMT_ASC_BASIC ) for ( int i=0; i<nTri; ++i) if ( tri[i].bdy ) printTriangle(&tri[i],level);
    
    if ( level >=FMT_ASC_BASIC ) cout <<" TETRAHEDRONS ON HULL\n=================\n";
    if ( level > FMT_ASC_BASIC ) for ( int i=0; i<nTtr; ++i) if ( ttr[i].bdy ) printTetrahedron(&ttr[i],level);
}

void Triangulation::printVertex(VertexType *v, int level) {
    
    cout << "Vrt " << v->idx << (v->bdy == true? 'H' : 'I') << "(";
    if (level>FMT_ASC_BASIC) for (int i2=0; i2<3; ++i2) cout << v->p[i2] <<";";
    cout <<  "):" << v->nVrtEdg << "Edg: " ;
    if (level>FMT_ASC_EXT) for (int i2=0; i2<v->nVrtEdg; ++i2) cout << v->edg[i2] <<";";
    cout << v->nVrtTri << "Tri: ";
    if (level>FMT_ASC_EXT)  for (int i2=0; i2<v->nVrtTri; ++i2) cout << v->tri[i2] <<";";
    cout << v->nVrtTtr << "Ttr: ";
    if ( level>FMT_ASC_EXT) for (int i2=0; i2<v->nVrtTtr; ++i2) cout << v->ttr[i2] <<";";
    cout << "\n" ;
}
void Triangulation::printEdge(EdgeType *e, int level) {
    cout << "Edg " << e->idx << (e->bdy == true? 'H' : 'I')<<"("<< e->vrt[0] << '-' <<  e->vrt[1] << ")"  <<":"<<  e->nEdgTri << "Tri: ";
    if ( level > FMT_ASC_EXT ) for    (int n=0; n<e->nEdgTri;   ++n) cout << e->tri[n] << " ";
    cout << ";  " << e->nEdgTtr << "Ttr: " ;
    if ( level > FMT_ASC_EXT ) for    (int n=0; n<e->nEdgTtr; ++n) cout << e->ttr[n] << " ";
    cout << "\n";
}

void Triangulation::printTriangle(TriType *t, int level) {
    
    cout <<"Tri " << t->idx << (t->bdy == true? 'H' : 'I') <<":(";
    if ( level >FMT_ASC_BASIC ) for ( int d=0; d<3; ++d ) cout << t->vrt[d] << ( d<2? "-" : ")"  );
    if ( level >FMT_ASC_EXT ) {
        cout << " Edg: ";
        for ( int d=0; d<3; ++d ) cout << t->edg[d] << ";";
        cout << " Ttr: ";
        for (int d=0; d<2; ++d) cout << t->ttr[d] << ";";
    }
    cout << '\n';
}

void Triangulation::printTetrahedron(TetraType *t, int level){
    cout << "Ttr " << t->idx << (t->bdy == true? 'H' : 'I') << "(";
    for (int i2=0; i2<3; ++i2) cout << t->c[i2] << " ";
    cout <<") Vrt:";
    for (int i2=0; i2<4; ++i2) cout << t->vrt[i2] << ";";
    if ( level > FMT_ASC_EXT) {
        cout << "Edg: ";
        if ( level>FMT_ASC_EXT ) for (int i2=0; i2<6; ++i2) cout << t->edg[i2] << ";";
        cout << "Tri: ";
        if ( level>FMT_ASC_EXT ) for (int i2=0; i2<4; ++i2) cout << t->tri[i2] << ";";
    }
    cout << "(" << t->n_adjacent << " Adj): ";
    if ( level>FMT_ASC_EXT ) for ( int i2=0;i2<t->n_adjacent; ++i2) cout << t->a[i2] << ";";
    cout << "(" << t->n_neighbor << " Ngh): ";
    if ( level>FMT_ASC_EXT ) for ( int i2=0;i2<t->n_neighbor; ++i2) cout << t->n[i2] << ";";
    cout << "\n";
}

void Triangulation::construct_normal(TetraType *t, TriType tri_loc,int ittr_loc){
    // construct outer normal of triangle with respect to tetrahedron t;
    VertexType *p_vrt0, *p_vrt1;
    double e_vec[3][3];
    
    for (int d=0; d<3; ++d )
    {
        p_vrt0 = &vrt[tri_loc.vrt[0]];
        p_vrt1 = &vrt[tri_loc.vrt[1]];
        e_vec[0][d]=p_vrt1->p[d] -p_vrt0->p[d];   // edge 0-1
        p_vrt1 = &vrt[tri_loc.vrt[2]];
        e_vec[1][d]=p_vrt1->p[d] -p_vrt0->p[d];   // edge 0-2
        e_vec[2][d]=tri_loc.c[d] - t->c[d];       // centroid triangle -- centroid tetrahedron
    }
    crs3(e_vec[0],e_vec[1],tri_loc.n[ittr_loc]);
    renorm3(tri_loc.n[ittr_loc]);
    renorm3(e_vec[2]);
    if ( dot3(tri_loc.n[ittr_loc],e_vec[2]) < 0. )  // turn around if inner normal
        for (int d=0;d<3;++d)
            tri_loc.n[ittr_loc][d]=-tri_loc.n[ittr_loc][d];
}


bool Triangulation::hasAtlas(){ return atlas.size() > 0? true : false; }

void Triangulation::setAtlas(string fname){
    int iLine=0, iAtl,iVrt, iTtr,reg_list[4],iuse;
    TetraType *pT;
    
    string line,dstr;
    ifstream fs;
    double dist[4];
    if ( file_exist(fname) )
        cout << "Setting Atlas from file " << fname << '\n';
    
    fs=ifstream(fname);
    // First, get region names

    atlas.push_back(fname);
    while ( getline(fs,line) ) {
        iLine++;
        vector<string> line_split;
        string_split(line, line_split,',');
        
        if ( iLine > nVrt ) {
            cout << "WARNING: Atlas of size " << atlas.size() << " contains line " << iLine << "; need " << nVrt << '\n';
            continue;
        }
        
        if ( find(atlas,line_split[1]) == atlas.size() )
            atlas.push_back(line_split[1]);
    }
    
    // Now, iterate over the file again and save label for each node
    fs.clear();  // clear ifstream in case, we have reached EOF
    fs.seekg(0,ios::beg);
    iLine=0;
    while ( getline(fs,line) ) {
        iLine++;
        
        if ( iLine > nVrt ) continue;
        
        vector<string> line_split;
        string_split(line,line_split,',');
        if ( (iAtl = find(atlas,line_split[1])) == atlas.size() ) {
            cout << "ERROR: Atlas type " << line_split[1] << "not found in atlas of size " << atlas.size() << '\n';
            exit(EXIT_FAILURE);
        } else {
            iVrt=atoi(line_split[0].c_str())-1;  // numbering in Atlas is 1-based;
            vrt[iVrt].atl=iAtl;
            vrt[iVrt].dia=(double)atof(line_split[3].c_str());
            vrt[iVrt].vol=(double)atof(line_split[2].c_str())/1000.;
        }
    }
    
    if (iLine < nVrt) {
        cout << "ERROR: Atlas of size " << atlas.size() << " only contains " << iLine << "nodes; need" << nVrt << '\n';
        // write msg
        // gbm_throw(msg,GBM_ATLAS_ERROR) -- write log / stdout / stack trace? / throw error
        exit(EXIT_FAILURE);
    }
    
    // Now set labels to tetrahedrons where all vertices belong to the same region
    for(iTtr=0; iTtr<nTtr; ++iTtr) {
        pT=&ttr[iTtr];
        pT->atl=NAN;

        for ( iVrt=0; iVrt<4; ++iVrt )
            reg_list[iVrt] = vrt[pT->vrt[iVrt]].atl;
        sort4(reg_list);
        
        if ( reg_list[0] == reg_list[3] || reg_list[0]==reg_list[2] ) // all or most ertices of the same region
            iuse=0;
        else if ( reg_list[1] == reg_list[3] )
            iuse=1;
        else
        {
            for (iVrt=0; iVrt<4; ++iVrt)
                dist[iVrt]=dist3(pT->c,vrt[pT->vrt[iVrt]].p);

            if( reg_list[0] == reg_list[1] && reg_list[2] == reg_list[3]) {// split in two
                if   ( dist[0]+dist[1] < dist[2]+dist[3] )
                    iuse=0;
                else
                    iuse=2;
            }
            else if ( reg_list[0] == reg_list[1] ||  reg_list[1] == reg_list[2] )
                iuse=1;
            else if ( reg_list[2] == reg_list[3] )
                iuse=2;
            else   // all four are different use region from closest vertex
                iuse=minpos(dist,4);
        }
        
        pT->atl=vrt[pT->vrt[iuse]].atl;
    }
    
    return;
}
