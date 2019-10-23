//  triangulate.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//
#include "triangulate.hpp"

Triangulation::Triangulation(){
/*    char fname[MAX_CHAR_LEN];
    
    cout << "Enter name of mesh file containing triangulation:";
    cin >> fname;
    
    read_mesh(fname);
*/
    return; 
}

Triangulation::Triangulation(char *fname,bool p[3]){
    nEdg=0;
    nTri=0;
    nTtr=0;
    nVrt=0;
    for (int id=0; id<3; ++id)
        periodic[id]=p[id];
    read_mesh(fname);
    
    return;
}

void Triangulation::read_mesh(char *fname){
    int iLine=0,nVrt_target;
    string line;
    istringstream lss;
    ifstream tfile;

    // check if file exists
    if ( file_exist(fname) )
        GBMLog("Reading mesh from file \'" + string( &fname[0] ) + "\'");
    
    tfile = ifstream(fname);
    
    // READ VERTICES
    getline(tfile,line); iLine++;
    lss=istringstream(line);
    lss >> nDim;
    if ( nDim == 4 )
        nDim--;   // decrement by one, because it includes elevation
    
    if ( nDim != 3 )
        GBMError("read_mesh", string("in"+ string(fname)+":"+ to_string(iLine+1)+
                                     "Dimensions of"+to_string(nDim) + " not supported"),GBMERR_DIM);
    getline(tfile,line); iLine++;
    lss = istringstream(line);  lss >> nVrt_target;
    GBMLog("Dimensions: " +to_string(nDim));
    GBMLog("Reading:    " +to_string(nVrt_target) +"Vertices");

    for ( int i=0; i<3; ++i ){
        box[2*i]= 9e9;
        box[2*i+1]=-9e9;
    }
    
    vrt.reserve(nVrt_target);
    for ( int iVrt=0; iVrt<nVrt_target; ++iVrt) {
        getline(tfile,line);
        if ( ( nVrt=addVertex(lss=istringstream(line)) ) != iVrt+1)
            GBMError("Triangulation::read_mesh",
                     "Reading Vertex "+to_string(iVrt)+" Position in File does not match position in Vector",
                     GBMERR_INDEX);
    }
    
    // READ TETRAHEDRONS
    getline(tfile,line); iLine++;
    lss=istringstream(line);
    lss >> nTtr;
    
    ttr.reserve(nTtr);
    volume = 0. ;
    GBMLog("Reading:    "+to_string(nTtr) +"Tetrahedrons");
    
    for ( int iTtr=0; iTtr<nTtr; ++iTtr) {
        getline(tfile,line);
        lss=istringstream(line);
        if ( addTetra(lss) != iTtr+1 )
            GBMError("Triangulation::read_mesh","Index of Tetrahedron " + to_string(iTtr) +
                            "does not match size of vector " +to_string(ttr.size()), GBMERR_INDEX);
    }
    
    fixNeighbors();
    GBMLog("Fixed Neighbors");
    
    // Construct Hull; fix boudary property of vertices, edges, and terahedras
    ConstructHull();
    GBMLog("Calculated Hull; contains " + to_string(nHul) + " triangles.");

    nTtr_inner=nTtr;
    nVrt_inner=nVrt;
    
    return;
}

void Triangulation::fixNeighbors() {
    TetraType *pt;
    int it, nEdgCommon, i_other, iAdj;
    
    for (it=0; it<nTtr; ++it) {
        pt=&ttr[it];
        pt->n_neighbor=0;        // reset number of neighbors

        // now check for all adjactent tetrahedrons whether a surface is shared
        for (iAdj=0; iAdj < pt->n_adjacent; iAdj++ ) {
            i_other = pt->a[iAdj];
            nEdgCommon=0;
            for ( int iEdg=0; iEdg<6; ++ iEdg)   // count all common edges
                if ( afind(ttr[i_other].edg,6,pt->edg[iEdg]) != 6 )
                    nEdgCommon++;
           
            if ( nEdgCommon == 3) {
                pt->n[pt->n_neighbor]=i_other;
                pt->nghVrt[pt->n_neighbor]=nonCommonVertex(it,i_other);
                pt->n_neighbor++;
                if ( pt->n_neighbor > 4 )
                    GBMError("Triangulation::fixNeighbor",
                             "Internal-the number of neighbors cannot be larger than 4 Found "+to_string(pt->n_neighbor) +"for Tetra " +to_string(it),
                             GBMERR_TOPOLOGY);
            } else if ( nEdgCommon != 1 )
                GBMError("Triangulation::fixNeighbor",
                         "Internal - the number of common Edges among Tetrahedrons Must either be 1 or 3, Found"
                         +to_string(nEdgCommon),GBMERR_TOPOLOGY);
        }
        if ( pt-> n_neighbor <1 || pt->n_neighbor > 4)
            GBMError("Triangulation::fixNeighbor",
                     "Internal -- the number of neighbors for a Tetrahedron must be 2,3, or 4; found "
                     +to_string(pt->n_neighbor), GBMERR_TOPOLOGY);
    }
    
    return;
}

int Triangulation::nonCommonVertex(int it1, int it2) {
    // returns the first vertex in ttr[it2] that is not contained in ttr[it1]
    int iv;
    int *t1=ttr[it1].vrt, *t2=ttr[it2].vrt;
    
    for ( iv=0; iv<4; ++iv )
        if ( t1[0]==t2[iv] or
             t1[1]==t2[iv] or
             t1[2]==t2[iv] or
             t1[3]==t2[iv])
            continue;
        else
            return t2[iv];
    return -1;
}

int Triangulation::addTris(int it, TetraType *t){
    int inew=0,v[3];
    int com[4][3] = { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}}; // 4 possible combinations of the three vertices
    int itri_loc;
    
    for ( int i=0; i<4; ++i){
        for (int dd=0; dd<3; ++dd) v[dd]=t->vrt[com[i][dd]];
        itri_loc=isTriangle(v);
        
        if ( itri_loc < 0) { // Triangle needs to be added
            TriType tri_loc;
            tri_loc.idx=nTri;
            itri_loc=nTri;
            tri_loc.edg[0]=isEdge(v[0],v[1]);
            tri_loc.edg[1]=isEdge(v[0],v[2]);
            tri_loc.edg[2]=isEdge(v[1],v[2]);
            
            for (int d=0; d<3; ++d) {
                tri_loc.vrt[d]=v[d];                 // add vertices
                tri_loc.c[d]=0.;                     // calculate centre
                for(int dd=0;dd<3;++dd)
                    tri_loc.c[d] += vrt[v[dd]].p[d] /3.;
                if ( vrt[v[d]].nVrtTri >= MAX_VRT_TRI )
                    GBMError("Triangulation::addTris",
                             "Too many triangles (" +to_string(vrt[v[d]].nVrtTri)+") for Vertex " +to_string(v[d]),
                             GBMERR_MAXVRTTRI);
                vrt[v[d]].tri[vrt[v[d]].nVrtTri]=itri_loc;
                vrt[v[d]].nVrtTri++;
            }
            
            tri_loc.ttr[0]=it;
            construct_normal(t,tri_loc,0);
            
            tri_loc.ttr[1]=-1;
            tri_loc.bdy=true;
            tri_loc.area=triArea(vrt[v[0]].p,vrt[v[1]].p,vrt[v[2]].p);
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
            int dum=(edg[tri[itri_loc].edg[d]].nEdgTri);
            if(  dum>= MAX_EDG_TRI ) {
                EdgeType *e=&edg[tri[itri_loc].edg[d]];
                GBMWarning("Triangulation::addTris",
                           "Edge "+to_string(tri[itri_loc].edg[d])+":"+
                           to_string(e->vrt[0]) + "-" + to_string(e->vrt[1]) +"("+
                           to_string(vrt[e->vrt[0]].p[0])+","+
                           to_string(vrt[e->vrt[0]].p[1])+","+
                           to_string(vrt[e->vrt[0]].p[2])+")-("+
                           to_string(vrt[e->vrt[1]].p[0])+","+
                           to_string(vrt[e->vrt[1]].p[1])+","+
                           to_string(vrt[e->vrt[1]].p[2])+")-(",
                           GBMERR_MAXEDGTRI);
                GBMError("Triangulation::addTris",string("Too many Triangles (" +to_string(dum) +") on Edge " +to_string(tri[itri_loc].edg[d])),
                         GBMERR_MAXEDGTRI);
            }
            edg[tri[itri_loc].edg[d]].tri[dum]=itri_loc;
            edg[tri[itri_loc].edg[d]].nEdgTri++;
        }
    }
    return inew;
}

int Triangulation::isTriangle(int v[3]){
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

int Triangulation::addEdges(int it, TetraType *t){
    int v0,v1, inew=0,iEdg,iEdg_loc=0,iTtr_loc;
    TetraType *t_other;
    VertexType *v;
    
    iEdg_loc=0;
    for ( int i0=0; i0<4; ++i0){
        for ( int i1=i0+1; i1<4; ++i1) {
            v0=t->vrt[i0];
            v1=t->vrt[i1];
            if ( v0 > v1 ) swap(v0,v1);
            iEdg = isEdge(v0,v1);
            if ( iEdg >= 0 )      // Edge exists already; add information about Tetrahedron
            {
                // recover other tetrahedrons from edge properties and tell them we are adjacent
                for ( int it_loc=0; it_loc<edg[iEdg].nEdgTtr; ++it_loc) {
                    iTtr_loc=edg[iEdg].ttr[it_loc];
                    if (afind(t->a,t->n_adjacent,iTtr_loc) == t->n_adjacent) {
                        if ( t->n_adjacent >= MAX_TTR_ADJ )
                            GBMError("Triangulation::addEdges",
                                     "Too many adjacent tetras (" +to_string(t->n_adjacent)+") on Tetra "
                                     +to_string(iTtr_loc),GBMERR_INDEX);
                        t->a[t->n_adjacent]=iTtr_loc;
                        t->n_adjacent++;
                    }
                    
                    t_other = &ttr[iTtr_loc];
                    
                    if ( afind(t_other->a,t_other->n_adjacent,it) == t_other->n_adjacent ) {
                        if ( t_other->n_adjacent >= MAX_TTR_ADJ )
                            GBMError("Triangulation::addEdges",
                                     "Too many adjacent tetras (" +to_string(t_other->n_adjacent)+") on Tetra" + to_string(iTtr_loc),GBMERR_INDEX);
                        t_other->a[t_other->n_adjacent]=it;
                        t_other->n_adjacent++;
                    }
                }
                if (edg[iEdg].nEdgTtr >= MAX_EDG_TTR )
                    GBMError("Triangulation::addEdges",
                             "Too many Tetras (" +to_string(edg[iEdg].nEdgTtr)+") on Edge "+to_string(iEdg),
                             GBMERR_MAXEDGTTR);
                             
                edg[iEdg].ttr[edg[iEdg].nEdgTtr]=it;
                edg[iEdg].nEdgTtr++;
                t->edg[iEdg_loc]=iEdg;
            }
            else                 // Edge does not exist yet; needs to be created
            {
                EdgeType edg_loc;
                edg.push_back(edg_loc);
                
                edg[nEdg].idx=nEdg;
                edg[nEdg].vrt[0]=v0;
                edg[nEdg].vrt[1]=v1;
                edg[nEdg].ttr[0]=it;
                edg[nEdg].nEdgTtr=1;
                edg[nEdg].nEdgTri=0;
                //edg.push_back(edg_loc);
                t->edg[iEdg_loc]=nEdg;
                
                v=&vrt[v0];
                if ( afind(v->edg,v->nVrtEdg,nEdg) == v->nVrtEdg) {
                    if ( v->nVrtEdg >= MAX_VRT_EDG )
                        GBMError("Triangulation::addEdges",
                                 "Too many Edges ("+to_string(v->nVrtEdg)+") on Vertex "+to_string(v0),
                                 GBMERR_MAXEDGVRT);
                    v->edg[v->nVrtEdg] = nEdg;
                    v->nVrtEdg++;
                }
                
                v=&vrt[v1];
                if ( afind(v->edg,v->nVrtEdg,nEdg) == v->nVrtEdg ){
                    if ( v->nVrtEdg >= MAX_VRT_EDG )
                        GBMError("Triangulation::addEdges",
                                 "Too many Edges (" +to_string(v->nVrtEdg)+") on Vertex "+to_string(v1),
                                 GBMERR_MAXEDGVRT);
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

int Triangulation::isEdge(int v0, int v1){
    // check if one of the vertices already holds this edge
    for ( int iEdg=0; iEdg < vrt[v0].nVrtEdg; ++iEdg ){
        if ( edg[vrt[v0].edg[iEdg]].vrt[1] == v1 )
            return vrt[v0].edg[iEdg];
    }
    return ( -1 );
}

int Triangulation::isVertex(double p[3], int start, int end){
    int iv,id;

    // we should use some pre-hashing, for instance according to
    // a global coarse cartesian grid here, to make the finding
    // faster. 
    
    
    for ( iv=start; iv<(end==-1? nVrt : end); ++iv) {
        for (id=0; id<3; ++id)
            if ( vrt[iv].p[id] != p[id])
                break;
        if ( id== 3 )
            return iv;
    }
 
    return -1;
}

int Triangulation::addVertex(istringstream &l, int halo){
    double p[4];
    l >> p[0] >> p[1] >> p[2] >> p[3];
    return addVertex(p,halo);
}
int Triangulation::addVertex(double p[4],int halo) {
    int i=int(vrt.size());
    double *b;
    VertexType v;
    
    vrt.push_back(v);
    for (int id=0; id<3; ++id) vrt[i].p[id]=p[id];
    
    if ( p[3] > 0 )
        vrt[i].e=p[3];
    else
        vrt[i].e=p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    
    vrt[i].idx=int(vrt.size()); // 0-based vertex indices
    vrt[i].bdy=false;
    vrt[i].nVrtEdg=0;
    vrt[i].nVrtTri=0;
    vrt[i].nVrtTtr=0;
    vrt[i].nVrtEdg=0;
    vrt[i].nVrtTri=0;
    vrt[i].nVrtTtr=0;
    vrt[i].halo=halo;
    
    if ( halo < 0 )
        b=&box[0];
    else
        b=&halo_box[0];
    
    for ( int id=0;id<3; ++ id) {
        b[2*id]   = vrt[i].p[id] < b[2*id]?   vrt[i].p[id]: b[2*id];    // min x
        b[2*id+1] = vrt[i].p[id] > b[2*id+1]? vrt[i].p[id]: b[2*id+1];  // max x
    }
    return int(vrt.size());
}

int Triangulation::isTetra(int p[4], int start, int end) {
    for(int it=start; it<(end==-1? nTtr : end); ++it)
        if (ttr[it].vrt[0]==p[0] and ttr[it].vrt[1]==p[1] and
            ttr[it].vrt[2]==p[2] and ttr[it].vrt[3]==p[3] )
            return it;
    return -1;
}

int Triangulation::addTetra(istringstream &l, int halo) {
    int p[4];
    l >> p[0] >> p[1] >> p[2] >> p[3];
    return addTetra(p,halo);
}

int Triangulation::addTetra(int p[4],int halo) {
    TetraType t;
    VertexType *pv;
    int iTtr = int(ttr.size());
    int id;

    for (id=0;id<4;++id) t.vrt[id]=p[id];
    sort4(t.vrt);
    t.idx=iTtr;
    t.n_adjacent=0;
    t.n_neighbor=0;
    t.bdy=false;
    t.halo=halo;
    
    // Calculate Centroid of tetrahedron
    for ( int idim=0; idim<3; ++idim) {
        t.c[idim]=0.;
        for ( int ivert=0; ivert<4; ++ivert) {
            t.c[idim] += vrt[t.vrt[ivert]].p[idim];
        }
        t.c[idim] *= 0.25;
    }
    
    for ( int i=0; i<4; ++i) {
        pv=&vrt[t.vrt[i]];
        if ( pv->nVrtTtr > MAX_VRT_TTR )
            GBMError("Triangulation::addTetra",
                     "Too many tetrahedrons ("+to_string(pv->nVrtTtr)+") on Vertex "+to_string(t.vrt[i]),
                     GBMERR_MAXTTRVRT); 
        pv->ttr[pv->nVrtTtr]=iTtr;
        pv->nVrtTtr++;
        t.n[i]=-1;
        t.nghVrt[i]=-1;
    }
    t.vol=tetraVolume(vrt[t.vrt[0]].p, vrt[t.vrt[1]].p,
                           vrt[t.vrt[2]].p, vrt[t.vrt[3]].p);
    volume += t.vol;
    addEdges(iTtr,&t);  // update list of edges
    addTris(iTtr,&t);   // update list of triangles
    
    ttr.push_back(t);

    return int(ttr.size());
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

    GBMLog("WRITING GRID: " + grid_file + "; Format_tag:" + to_string(format));
    
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
        if ( hasAtlas() ) {
        // <DataArray type="Int32" Name="VrtAtlasRegion" format="ascii">
        att.push_back("type");  att.push_back("Int32");
        att.push_back("Name");  att.push_back("VrtAtlasRegion");
        att.push_back("format");att.push_back("ascii");
        for(iv=0; iv<nVrt;++iv )
        {
            idata[iv]    =vrt[iv].atl;
            data[iv]     =vrt[iv].dia;
            data[nVrt+iv]=vrt[iv].vol;
        }
        vtkXMLWriteDataArray(gfile,att,nVrt,idata);
        att[1]="Float32"; att[3]="VrtAtlasDiameter";
        vtkXMLWriteDataArray(gfile,att,nVrt,data);
        
        att[3]="VrtAtlasVolume";
        vtkXMLWriteDataArray(gfile, att, nVrt, &data[nVrt]);
        }
        att.clear();
        
        att.push_back("type");  att.push_back("Int32");
        att.push_back("Name");  att.push_back("nVrtTri");
        att.push_back("format");att.push_back("ascii");
        for ( iv=0;iv<nVrt; ++iv )
            idata[iv] = vrt[iv].nVrtTri;
        vtkXMLWriteDataArray(gfile,att,nVrt,idata);
       
        att[3]="nVrtEdg";
        for ( iv=0;iv<nVrt; ++iv )
            idata[iv] = vrt[iv].nVrtEdg;
        vtkXMLWriteDataArray(gfile,att,nVrt,idata);
        
        att[3]="nVrtTtr";
        for ( iv=0;iv<nVrt; ++iv )
            idata[iv] = vrt[iv].nVrtTtr;
        vtkXMLWriteDataArray(gfile,att,nVrt,idata);
        
        att[3]="VrtHalo";
        for ( iv=0;iv<nVrt; ++iv )
            idata[iv] = vrt[iv].halo;
        vtkXMLWriteDataArray(gfile,att,nVrt,idata);
        
        gfile << "</PointData>\n";
        
        gfile << "<CellData>\n";
        // <DataArray type="Int32" Name="Atlas" format="ascii">
        att[1]="Int32"; att[3]="TtrAtlas";
        for ( it=0; it<nTtr; ++it) {
            idata[it]=ttr[it].atl;
            data[it] =ttr[it].vol;
        }
        vtkXMLWriteDataArray(gfile, att, nTtr, idata);
        
        att[1]="Float32";att[3]="TtrVolume";
        vtkXMLWriteDataArray(gfile, att, nTtr, data);
            
        att[1]="Int32"; att[3]="TtrHalo";
        for (it=0; it<nTtr; ++it)
            idata[it]=ttr[it].halo;
        vtkXMLWriteDataArray(gfile, att, nTtr, idata);
       
        att[3]="nTtrN_neighbor";
        for (it=0; it<nTtr; ++it)
            idata[it]=ttr[it].n_neighbor;
        vtkXMLWriteDataArray(gfile, att, nTtr, idata);
            
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
        GBMLog("Setting Atlas from file "+fname);
    
    fs=ifstream(fname);
    // First, get region names

    atlas.push_back(fname);
    while ( getline(fs,line) ) {
        iLine++;
        vector<string> line_split;
        string_split(line, line_split,',');
        
        if ( iLine > nVrt ) {
            GBMWarning("Triangulation::setAtlas",
                       "Atlas of size " + to_string(atlas.size()) + " contains line " +to_string(iLine)+"; need "+to_string(nVrt),
                       GBMERR_INDEX);
            continue;
        }
        
        if ( find(atlas,line_split[1]) == atlas.size() )
            atlas.push_back(line_split[1]);
    }
    
    // Now, iterate over the file again and save an integer label for each node
    fs.clear();           // clear ifstream in case, we have reached EOF
    fs.seekg(0,ios::beg); // go back to beginning of file
    iLine=0;
    while ( getline(fs,line) ) {
        iLine++;
        
        if ( iLine > nVrt ) continue;
        
        vector<string> line_split;
        string_split(line,line_split,',');
        if ( (iAtl = find(atlas,line_split[1])) == atlas.size() )
            GBMError("Triangulation::setAtlas",
                     "Atlas type " + line_split[1] +" not found in atlas of size " + to_string(atlas.size()),
                     GBMERR_INDEX);
        else {
            iVrt=atoi(line_split[0].c_str())-1;  // numbering in Atlas is 1-based;
            vrt[iVrt].atl=iAtl;
            vrt[iVrt].dia=(double)atof(line_split[3].c_str());
            vrt[iVrt].vol=(double)atof(line_split[2].c_str())/1000.;
        }
    }
    
    if (iLine < nVrt)
        GBMError("Triangulation::setAtlas",
                 "Atlas of size " + to_string(atlas.size()) + " only contains " +to_string(iLine)+"nodes; need" +to_string(nVrt),
                 GBMERR_INDEX);
    
    // Now estimate labels for tetrahedra from vertices
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

int Triangulation::TtrSplinesND(int sType) {
    int nd=-1;
    switch ( sType ) {
        case SPLINE_O1:
            nd=4;  break;
        case SPLINE_O2_LOCAL:
        case SPLINE_O2_EDGES:
        case SPLINE_O2_CLSTR:
            nd=10; break;
        default:
            GBMError("Triangulation::TtrSplinesAlloc",
                     "Unknown Spline Type " + to_string(sType),
                     GBMERR_PARAMETER);
    }
    return nd;
}

void Triangulation::TtrSplinesAlloc(int sType){
    int id,id2,iptr;
    int QRsize_mat, QRsize_vec, QRsize_bare;
    const int nd=TtrSplinesND(sType);
    
    QRsize_mat = nTtr*nd*nd;
    QRsize_vec = nTtr*nd;
    QRsize_bare= QRsize_mat + 3*QRsize_vec;
    //
    QRSplines_bare[sType] = (double *)malloc(QRsize_bare*sizeof(double));
    QRSplines_a[sType]= (double ***) malloc(nd*sizeof(double **));
    QRSplines_b[sType] = (double **)  malloc(nd*sizeof(double *));
    QRSplines_c[sType] = (double **)  malloc(nd*sizeof(double *));
    QRSplines_d[sType] = (double **)  malloc(nd*sizeof(double *));
    
    iptr=0;
    for ( id=0; id<nd; ++id ) {
        QRSplines_a[sType][id]= (double **) malloc(nd*sizeof(double*));
        for ( id2=0;id2<nd; ++id2) {
            QRSplines_a[sType][id][id2]=&(QRSplines_bare[sType][iptr]); iptr += nTtr;
        }
    }
    
    for ( id=0; id<nd; ++id ) {
        QRSplines_b[sType][id] = &(QRSplines_bare[sType][iptr]);  iptr += nTtr;}
    for ( id=0; id<nd; ++id ) {
        QRSplines_c[sType][id] = &(QRSplines_bare[sType][iptr]);  iptr += nTtr;}
    for ( id=0; id<nd; ++id ) {
        QRSplines_d[sType][id] = &(QRSplines_bare[sType][iptr]);  iptr += nTtr;}
    GBMLog("Prepared Left-Hand Side for Spline Calculations (nd="+to_string(nd)+")");

}

void Triangulation::TtrSplinesLHS(int sType){
    int iTtr, id;
    double dp[3];
    const int nd=TtrSplinesND(sType);
    int pts_O2_LOCAL[10][2]= SPLINE_O2_LOCAL_POINTS;
    int pts_O2_EDGES[10][2]= SPLINE_O2_LOCAL_POINTS;
    int pts[10][2];
    // calculate spline
    switch ( sType ) {
        case SPLINE_O1: {
            VertexType *v;
            for (id=0;id<nd; ++id) {
                for (iTtr=0; iTtr<nTtr; ++iTtr) {
                    v=&(vrt[ttr[iTtr].vrt[id]]);
                    QRSplines_a[sType][0][id][iTtr]=1.;
                    QRSplines_a[sType][1][id][iTtr]=v->p[0] - ttr[iTtr].c[0];
                    QRSplines_a[sType][2][id][iTtr]=v->p[1] - ttr[iTtr].c[1];
                    QRSplines_a[sType][3][id][iTtr]=v->p[2] - ttr[iTtr].c[2];
                }
            }
            break; }
        case SPLINE_O2_LOCAL: {
            memcpy(pts,pts_O2_LOCAL,20*sizeof(int)); }
        case SPLINE_O2_EDGES: {
            if (sType==SPLINE_O2_EDGES) memcpy(pts,pts_O2_EDGES,20*sizeof(int));
            VertexType *v[2];
            for(id=0;id<nd; ++id) {
                for (iTtr=0; iTtr<nTtr_inner; ++iTtr) {
                    if ( pts[id][0] < 4 )
                        v[0]=&vrt[ttr[iTtr].vrt[pts[id][0]]];
                    else
                        v[0]=&vrt[ttr[iTtr].nghVrt[pts[id][0]-4]];
                    if ( pts[id][1] < 4 )
                        v[1]=&vrt[ttr[iTtr].vrt[pts[id][1]]];
                    else
                        v[1]=&vrt[ttr[iTtr].nghVrt[pts[id][1]-4]];
                    
                    QRSplines_a[sType][0][id][iTtr]=1.;
                    for ( int idim=0; idim<3; ++idim) {
                        dp[idim]=(v[0]->p[idim]+v[1]->p[idim])/2. - ttr[iTtr].c[idim];
                        QRSplines_a[sType][1+idim][id][iTtr]=dp[idim];
                        QRSplines_a[sType][4+idim][id][iTtr]=(dp[idim]*dp[idim])/2.;
                    }
                    QRSplines_a[sType][7][id][iTtr]=dp[0]*dp[1]/2.;
                    QRSplines_a[sType][8][id][iTtr]=dp[0]*dp[2]/2.;
                    QRSplines_a[sType][9][id][iTtr]=dp[1]*dp[2]/2.;
                }
            }
            break; }
        default:
            GBMError("Triangulation::TtrSplinesLHS",
                     "SplineType " + to_string(sType) + " not implemented",
                     GBMERR_UNDEVELOPED);
    }
    m_in_qrdcmp(QRSplines_a[sType], nd, nd, nTtr, QRSplines_c[sType], QRSplines_d[sType]);
    return;
}

void Triangulation::TtrSplinesRHS(int sType, double *data){
    // data contains vertex data;
    int iTtr,id;
    const int nd=TtrSplinesND(sType);

    switch ( sType ) {
        case SPLINE_O1: {
            for(id=0;id<4;++id)
                for(iTtr=0;iTtr<nTtr;++iTtr)
                    QRSplines_b[sType][id][iTtr] = data[ttr[iTtr].vrt[id]];
            break; }
        case SPLINE_O2_LOCAL: {
            int v[2];
            int pts[10][2] = SPLINE_O2_LOCAL_POINTS;
            for(id=0;id<nd;++id)
                for ( iTtr=0;iTtr<nTtr; ++iTtr) {
                    v[0]=ttr[iTtr].vrt[pts[id][0]];
                    v[1]=ttr[iTtr].vrt[pts[id][1]];
                    QRSplines_b[sType][id][iTtr] = (data[v[0]]+data[v[1]])/2. ;
                }
            break; }
        case SPLINE_O2_EDGES: {
            int v[2];
            int pts[10][2] = SPLINE_O2_EDGES_POINTS;
            for(id=0;id<nd;++id)
                for ( iTtr=0;iTtr<nTtr_inner;++iTtr) {
                    if ( pts[id][0] < 4 )
                        v[0]=ttr[iTtr].vrt[pts[id][0]];
                    else
                        v[0]=ttr[iTtr].nghVrt[pts[id][0]-4];
                    if ( pts[id][1] < 4 )
                        v[1]=ttr[iTtr].vrt[pts[id][1]];
                    else
                        v[1]=ttr[iTtr].nghVrt[pts[id][1]-4];
                    if ( v[0] < 0 or v[1] <0 )
                        GBMError("Triangulation::TtrSplinesRHS","Neighbor Vertex does not exist",GBMERR_HALO); 
                    QRSplines_b[sType][id][iTtr] = (data[v[0]]+data[v[1]])/2. ;
                }
            break; }
        default:
            GBMError("Triangulation::TtrSplinesRHS",
            "SplineType " + to_string(sType) + " not implemented",
            GBMERR_UNDEVELOPED);
    }
    
    m_in_qrsolv(QRSplines_a[sType], nd, nd, nTtr,
                QRSplines_c[sType], QRSplines_d[sType], QRSplines_b[sType]);
    
    return;
}

void Triangulation::TtrCentroidSplines(int sType, double *v)
{
    // Evaluate Splines at Centroids of tetrahedrons
    // all splines have the centroid value as first entry
    double max=-9.e30;
    double min= 9.e30;
    for (int iTtr=0; iTtr<nTtr; ++iTtr) {
        v[iTtr] = QRSplines_b[sType][0][iTtr];
        if ( v[iTtr] > max ) max=v[iTtr];
        if ( v[iTtr] < min ) min=v[iTtr];
    }
    //cout << sType << "Centroids ("<<sType << "): [" << min << " , " << max << "]." << std::endl;
    return;


}


void Triangulation::TtrDerivativeSplines(int sType, int direction,double *v){
    
    double min, max;
    max=-9.e30;
    min= 9.e30;
    
    for (int iTtr=0; iTtr<nTtr; ++iTtr) {
        v[iTtr] = QRSplines_b[sType][direction][iTtr];
        if ( v[iTtr] > max ) max=v[iTtr];
        if ( v[iTtr] < min ) min=v[iTtr];
    }
    //cout << sType << "Derivative " << sType << "Dir" << direction << ": [" << min << " , " << max << "]." << std::endl;
    return;
}


void Triangulation::ConstructHull(){
    TriType *pt;
    for ( int iTri=0; iTri<nTri; ++iTri)
       if ( tri[iTri].bdy == true )
       {
           pt=&tri[iTri];
           hul.push_back( pt );
           for (int d=0; d<3; ++d) {
               edg[pt->edg[d]].bdy=true;
               vrt[pt->vrt[d]].bdy=true;
           }
           ttr[pt->ttr[0]].bdy=true;
           if ( pt->ttr[1] > 0 )
               GBMError("Triangulation::ConstructHull",
                        "Inconsistent boundary property on Triangle" + to_string(iTri)
                        + "(" +to_string(pt->vrt[0]) + "-" + to_string(pt->vrt[1]) +"-"+to_string(pt->vrt[2]),
                        GBMERR_HALO);
           nHul++;
       }
    return;
}

void Triangulation::ConstructHalo(){
    double x[4],x2[3],s[3],tmp1[4];
    vector<VertexType> hVrt;
    vector<EdgeType>   hEdg;
    vector<TetraType>  hTtr;
    vector<TriType>    hTri;
    vector<int>        vList;
    vector<bool>       vper[3];
    vector<vector<int>>newTtr;
    vector<vector<int>>newVrt;
    int iv,iv2,id,id2,ie,it,v_add,ttr_new[5];
    int nBdy=0;
    int dimCount=0;
    bool per[3],per_loc[3];
    
    GBMLog("Constructing Halo Region");
    
    nVrt_inner=nVrt;
    nTtr_inner=nTtr;

    tmp1[3]=-1.;
    for (id=0;id<3;++id)
        s[id] = fabs(box[2*id+1]-box[2*id]);
    
    for ( iv=0; iv<nVrt_inner; ++iv) {
        if ( vrt[iv].bdy ) {
            
            for (id=0;id<3;++id)
                   x[id]=vrt[iv].p[id];

            for ( iv2=0; iv2<nVrt_inner; ++iv2 ) {
                for (id=0; id<3; ++id)
                    x2[id] = vrt[iv2].p[id];

                for ( id=0; id<3; ++id){
                    if (( x2[id] == x[id]+s[id] or x2[id] == x[id]-s[id] )  and
                        ( x2[0]-x[0]==0 or fabs(x2[0]-x[0])==s[0]) and
                        ( x2[1]-x[1]==0 or fabs(x2[1]-x[1])==s[1]) and
                        ( x2[2]-x[2]==0 or fabs(x2[2]-x[2])==s[2]) ) {
                        // Any point on boundary (includes edges and corners)
                        for ( ie=0; ie<vrt[iv2].nVrtEdg; ++ie) { // loop over edges corresponding to p
                            if ( (v_add=edg[vrt[iv2].edg[ie]].vrt[0]) == iv2 )
                                v_add=edg[vrt[iv2].edg[ie]].vrt[1];
                            for ( id2=0;id2<3;++id2) tmp1[id2]=vrt[v_add].p[id2];
                            tmp1[id] += tmp1[id]>=x[id]? -s[id] : s[id];
                            if ( tmp1[id] < box[2*id] or tmp1[id] > box[2*id+1] )
                                if ( isVertex(tmp1,nVrt_inner)<0) nVrt=addVertex(tmp1,v_add);
                        }
                        
                        for ( it=0; it<vrt[iv2].nVrtTtr; ++it ){
                            ttrPeriodic(it,iv2,id,ttr_new);
                            newTtr.push_back(vector<int>(ttr_new,ttr_new+5));
                        }
                        
                        // Now treat the edge points (2 periodic dimensions)
                        dimCount=0;
                        for ( id2=0; id2<3; ++id2)
                            if ( fabs(x2[id2]-x[id2])==s[id2] and dimCount < 2 ) {
                                per[id2]=true;
                                dimCount++;
                            } else per[id2]=false;
                        
                        if ( dimCount == 2 ) {
                            for ( ie=0; ie<vrt[iv2].nVrtEdg; ++ie) {
                                for (id2=0;id2<3;++id2) per_loc[id2]=false;
                                if ( ( v_add=edg[vrt[iv2].edg[ie]].vrt[0]) == iv2 )
                                    v_add=edg[vrt[iv2].edg[ie]].vrt[1];
                                dimCount=0;
                                for ( id2=0;id2<3;++id2) {
                                    tmp1[id2]=vrt[v_add].p[id2];
                                    if ( per[id2] and dimCount <2 ) {
                                        tmp1[id2] += tmp1[id2]>=x[id2]? -s[id2] : s[id2];
                                        per_loc[id2]=true;
                                        ++dimCount;
                                    }
                                }
                                if ( isVertex(tmp1,0) < 0 ) nVrt=addVertex(tmp1,v_add);
                            }
                            for ( it=0; it<vrt[iv2].nVrtTtr; ++it ){
                                ttrPeriodic(it,iv2,per,ttr_new);
                                newTtr.push_back(vector<int>(ttr_new,ttr_new+5));
                            }
                        }
                    }
                } // for ( id=0; id<3; ++id)
                
                // Now treat the corner points (3 periodic dimensions)
                for ( id2=0; id2<3; ++id2)
                    if ( fabs(x2[id2]-x[id2])==s[id2] ) per[id2]=true;
                    else                                per[id2]=false;
                        
                if ( per[0] and per[1] and per[2] ) {
                    for ( ie=0; ie<vrt[iv2].nVrtEdg; ++ie) {
                        if ( ( v_add=edg[vrt[iv2].edg[ie]].vrt[0]) == iv2 )
                            v_add=edg[vrt[iv2].edg[ie]].vrt[1];
                        for ( id2=0;id2<3;++id2) {
                            tmp1[id2]=vrt[v_add].p[id2];
                            tmp1[id2] += tmp1[id2]>=x[id2]? -s[id2] : s[id2];
                        }
                        if ( isVertex(tmp1,0) < 0 ) nVrt=addVertex(tmp1,v_add);
                    }
                    for ( it=0; it<vrt[iv2].nVrtTtr; ++it ){
                        ttrPeriodic(it,iv2,per,ttr_new);
                        newTtr.push_back(vector<int>(ttr_new,ttr_new+5));
                    }
                }
            } // for (iv2=0; iv2<nVrt_inner; ++iv2)
            ++nBdy;
        } // if (vrt[iv].bdy)
    } // for (iv=0; iv<nVrt_inner; ++iv)
 
    // Finally, add the Tetrahedrons for the halo
    {   // performance-critical part, therefore instrumented.
        uint64_t ttime1=0,ttime2=0,ctime=0;
        int iloc;
        GBMLog("Adding " +to_string(newTtr.size()) + " Tetrahedra for Halo region");
        for ( it=0; it<newTtr.size(); ++it){
            ctime=current_time();
            // we know, it is on the halo, so we can skip the inner part when searching
            iloc = isTetra(newTtr[it].data(),nTtr_inner);
            ttime1 += current_time()-ctime;
            if ( iloc  < 0 ) {
                ctime=current_time();
                nTtr=addTetra(newTtr[it].data(),newTtr[it][4]);
                ttime2+= current_time()-ctime;
            }
        }
        GBMLog("...timing for isTetra: " +to_string(ttime1/1.e3) +" ms");
        GBMLog("...timing for addTetra:" +to_string(ttime2/1.e3) +" ms");
    }
    // Recalculate neighbor topology --> makes halo known to other nodes
    fixNeighbors();
    
    GBMLog("...Constructed Halo containing " +to_string(nVrt-nVrt_inner) + "vertices." );
    
    return;
}


void Triangulation::ttrPeriodic(int it, int iv2, int id, int *ttr_new){
    bool per[3] = {false, false, false};
    per[id]=true;
    ttrPeriodic(it,iv2,per, ttr_new);
}

void Triangulation::ttrPeriodic(int it, int iv2, bool per[3], int *ttr_new){
    int iv,id2, ttr_old[4];
    double min,s,tmp[4][3];
    bool beg;
     
    ttr_new[4]=it;
    
    for(iv=0;iv<4;++iv) {
        for (id2=0;id2<3;++id2)
            tmp[iv][id2]=vrt[ttr[vrt[iv2].ttr[it]].vrt[iv]].p[id2];
        ttr_old[iv]=isVertex(tmp[iv]);
    }
      
    for (id2=0;id2<3;++id2)
        if ( per[id2] ) {
            s=box[2*id2+1]-box[2*id2];
            min=9e9;
            for(iv=0;iv<4;++iv)
                min=tmp[iv][id2]<min? tmp[iv][id2] : min;
            beg = min < s/2.? true : false;
            for(iv=0;iv<4;++iv)
                tmp[iv][id2] += beg? s: -s;
        }

    // find first vertex
    if (( tmp[0][0] < box[0] or tmp[0][0] > box[1] ) or
        ( tmp[0][1] < box[2] or tmp[0][1] > box[3] ) or
        ( tmp[0][2] < box[4] or tmp[0][2] > box[5] ) )
        ttr_new[0]=isVertex(tmp[0],nVrt_inner);  // point lies on halo
    else
        ttr_new[0]=isVertex(tmp[0],0,nVrt_inner);// point in interior
    
    for ( iv=1;iv<4; ++iv) {
        int v_ref, ie, ne;
        bool found=false;
        
        ttr_new[iv]=-1;
        
        // check if we can identify further vertices as neighbors of the previous ones
        ne=vrt[ttr_new[0]].nVrtEdg;
        for ( iv2=0; iv2< iv; ++iv2) {
            for ( ie=0; ie<ne; ++ie) {
                if ( (v_ref=edg[ vrt[ttr_new[0]].edg[ie] ].vrt[0]) ==  ttr_new[0] )
                    v_ref=edg[ vrt[ttr_new[0]].edg[ie] ].vrt[1];
                if (vrt[v_ref].p[0] == tmp[iv][0] and
                    vrt[v_ref].p[1] == tmp[iv][1] and
                    vrt[v_ref].p[2] == tmp[iv][2] ) {
                    ttr_new[iv]=v_ref;
                    found=true;
                    break;
                }
            }
            if ( ie < ne ) break;
        }
        
        if ( ! found ) {
            if (( tmp[iv][0] < box[0] or tmp[iv][0] > box[1] ) or
                ( tmp[iv][1] < box[2] or tmp[iv][1] > box[3] ) or
                ( tmp[iv][2] < box[4] or tmp[iv][2] > box[5] ) )   // point lies on halo; skip interior
                ttr_new[iv]=isVertex(tmp[iv],nVrt_inner);
            else                                                   // point in interior
                ttr_new[iv]=isVertex(tmp[iv]);
        }
    
        if ( ttr_new[iv] < 0 )
            GBMError("Triangulation::ttrPeriodic","Vertex for Periodic Tetrahedron not found",GBMERR_HALO);
    }
    sort4(ttr_new);
    return ;
}
