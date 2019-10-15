//
//  vspline.cpp
//  vspline
//
//  Created by Cedrick Ansorge on 10.09.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//
#include "vspline.hpp"
int main(int argc, const char * argv[]) {
    
    double *dat_v, *dat_t;
    double k[3],l2,linf,l2_ref;
    string ofile, gridFormat;
    string testGrid,rundir;
    int testList[] = {2,4,6,8,10,12,16};
    int nTest = 7;

    Namelist nl("/Users/zrick/WORK/research_projects/GBM/vspline.nml");
    nl.getList_dbl("Valid","k_test",&k[0],3);
    rundir=nl.getVal_str("Main","rundir");
    ofile =nl.getVal_str("Valid","output");
    
    for (int iTest=0; iTest<nTest; ++iTest) {
        string testGrid("/Users/zrick/Work/research_projects/GBM/cgal_3periodic/triangulation_");
        testGrid += std::to_string(testList[iTest]);
        GBMLog(string(to_string(iTest) + ":" +testGrid));
        
        Triangulation tri((char *)&testGrid.c_str()[0]);
        
        dat_v = (double *)malloc(tri.nVrt_inner*sizeof(double));
        dat_t = (double *)malloc(tri.nTtr_inner*sizeof(double));
    
        for (int i=0;i<tri.nVrt_inner;++i)
            dat_v[i]=test_func(tri.vrt[i].p,&k[0]);
    
        for (int i=0;i<tri.nTtr_inner;++i)
            dat_t[i]=0.;
    
        tri.TtrSplinesAlloc();
        tri.TtrSplinesLHS();
        tri.TtrSplinesRHS(dat_v);
        tri.CentroidSplines(dat_t);
    
        l2=0;
        linf=0;
        for (int iTtr=0; iTtr<tri.nTtr_inner; ++iTtr){
            double tmp,tmp2;
            tmp2= test_func(tri.ttr[iTtr].c,k);
            tmp = fabs(dat_t[iTtr] - tmp2);
            l2 += tmp*tmp;
            l2_ref += tmp2*tmp2;
            linf = ( tmp>linf ) ? tmp : linf;
        }
        cout << testList[iTest] << "," << l2/tri.nTtr << "," << l2/l2_ref << "," << linf << std::endl;
        writeSpline_test(ofile,&tri, dat_v, dat_t);
        
        free(dat_v);
        free(dat_t);
    }
    
    exit(EXIT_SUCCESS);
    
}

double test_func(double *p,double *k){
    return (cos(2*PI*k[0]*p[0])*cos(2*PI*k[1]*p[1])*cos(2*PI*k[2]*p[2]));
}


void writeSpline_test(string ofile, Triangulation *tri, double *dat_v, double *dat_t) {

    double *data;
    int *idata;
    int iv,it,d,data_size;
    int nTtr, nVrt;
    vector<string> att;
    ofstream gfile;
    
    nTtr=tri->nTtr;
    nVrt=tri->nVrt;
    
    data_size = 3*nVrt > 4*nTtr? 3*nVrt : 4*nTtr;
    
    data =(double*)malloc(data_size*sizeof(double));
    idata=(int*)   malloc(data_size*sizeof(int));
    
    vtkXMLFileOpen(ofile,nVrt,nTtr,gfile);
    
    // GRID INFORMATION -- VERTICES
    gfile << "<Points>\n";
    // <DataArray type="Float32" NumberOfComponents="3" format="ascii">
    att.push_back("type");              att.push_back("Float32");
    att.push_back("NumberOfComponents");att.push_back("3");
    att.push_back("format");            att.push_back("ascii");
    
    for (iv=0;iv<nVrt;++iv)
        for (d=0;d<3;++d)
            data[3*iv+d]=tri->vrt[iv].p[d];
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
            idata[4*it+d] = tri->ttr[it].vrt[d];
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
    
    // <DataArray type="Float32" Name="value_vrt" format="ascii">
    att.push_back("type");  att.push_back("Float32");
    att.push_back("Name");  att.push_back("value_vrt");
    att.push_back("format");att.push_back("ascii");

    vtkXMLWriteDataArray(gfile,att,nVrt,dat_v);

    gfile << "</PointData>\n";
    
    gfile << "<CellData>\n";
    // <DataArray type="Float32" Name="value_ttr" format="ascii">
    att[3]="value_ttr";
    vtkXMLWriteDataArray(gfile, att, nTtr, dat_t);
    gfile << "</CellData>\n";
    
    // FOOTER of VTK XML FILE
    vtkXMLFileClose(gfile);
    
    att.clear();
    free(idata);
    free(data);
}
