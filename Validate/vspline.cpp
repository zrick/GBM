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
    double k[3],l2,linf,l2_ref,l2_cube, linf_cube, l2_ref_cube,vref,vfil;
    string ofile, gridFormat;
    string testGrid,rundir;
    int testList[] = {2,3,4,6,8,10,12,16,20,30};
    int nTest = 10;//10;
    int sType = SPLINE_O1;
    int ii,jj,kk,nn;
    double x0[3], ***val3d;
    bool p[3]={false,false,false};

    Namelist nl("/Users/zrick/WORK/research_projects/GBM/vspline.nml");
    nl.getList_dbl("Valid","k_test",&k[0],3);
    rundir=nl.getVal_str("Main","rundir");
    ofile =nl.getVal_str("Valid","output");
    
    for (int iTest=0; iTest<nTest; ++iTest) {
        string testGrid("/Users/zrick/Work/research_projects/GBM/cgal_3periodic/triangulation_");
        testGrid += std::to_string(testList[iTest]);
        GBMLog(string(to_string(iTest) + ":" +testGrid));
        
        Triangulation tri((char *)&testGrid.c_str()[0],p);
        p[0]=true;
        p[1]=true;
        p[2]=true; 
        tri.ConstructHalo();
        
        dat_v = (double *)malloc(tri.nVrt*sizeof(double));
        dat_t = (double *)malloc(tri.nTtr*sizeof(double));
    
        for (int i=0;i<tri.nVrt;++i)
            dat_v[i]=test_func(tri.vrt[i].p,&k[0]);
    
        for (int i=0;i<tri.nTtr;++i)
            dat_t[i]=0.;
    
        tri.TtrSplinesAlloc(sType);
        tri.TtrSplinesLHS(sType);
        tri.TtrSplinesRHS(sType,dat_v);
        tri.TtrCentroidSplines(sType,dat_t);
    
        l2=0.;
        l2_ref=0.;
        linf=0.;
        
        for (int iTtr=0; iTtr<tri.nTtr_inner; ++iTtr){
            double tmp,tmp2;
            tmp2= test_func(tri.ttr[iTtr].c,k);
            tmp = fabs(dat_t[iTtr] - tmp2);
            l2 += tmp*tmp;
            l2_ref += tmp2*tmp2;
            linf = ( tmp>linf ) ? tmp : linf;
        }
        
        nn=testList[iTest];
        val3d= (double ***) malloc ( (nn+1) * sizeof(double**) );
        for (ii=0; ii<=nn; ++ii) {
            x0[0]=ii/double(nn);
            val3d[ii] = (double **) malloc ( (nn+1)*sizeof(double *));
            for (jj=0;jj<=nn; ++jj) {
                val3d[ii][jj] = (double *) malloc( (nn+1)*sizeof(double));
                x0[1]=jj/double(nn);
                for ( kk=0; kk<=nn; ++kk) {
                    x0[2] = kk/double(nn);
                    val3d[ii][jj][kk] = test_func(x0,k);
                }
            }
        }

        l2_cube=0.;
        l2_ref_cube=0.;
        linf_cube=0.;
        for (ii=0; ii<nn; ++ii) {
            x0[0]=(ii+0.5)/double(nn);
            for (jj=0;jj<nn; ++jj) {
                x0[1]=(jj+0.5)/double(nn);
                for ( kk=0; kk<nn; ++kk) {
                    x0[2] = (kk+0.5)/double(nn);
                    vref=test_func(x0,k);
                    vfil=(val3d[ii  ][jj  ][kk] + val3d[ii  ][jj  ][kk+1]
                        +val3d[ii  ][jj+1][kk] + val3d[ii  ][jj+1][kk+1]
                        +val3d[ii+1][jj  ][kk] + val3d[ii+1][jj  ][kk+1]
                        +val3d[ii+1][jj+1][kk] + val3d[ii+1][jj+1][kk+1])/8.;
                    l2_cube += (vfil-vref)*(vfil-vref);
                    linf_cube = fabs(vfil-vref) > linf_cube? fabs(vfil-vref) : linf_cube;
                    l2_ref_cube += vref*vref;
                }
            }
        }
        
        GBMLog(to_string(testList[iTest])+","+to_string(l2/tri.nTtr)+","+to_string(l2/l2_ref)+","+to_string(linf)+
               "|"+to_string(l2_cube/(nn*nn*nn))+","+to_string(l2_cube/l2_ref_cube)+","+to_string(linf_cube));
        cout << testList[iTest] << " " << l2/tri.nTtr << " " <<  l2/l2_ref << " " << linf << "|";
        cout << testList[iTest] << " " << l2_cube/(nn*nn*nn) << " " <<  l2_cube/l2_ref_cube << " " << linf_cube << std::endl;
        stringstream sstream;
        sstream << ofile << "_" << testList[iTest] << ".vtu.xml";
        writeSpline_test(sType, sstream.str(),&tri, dat_v, dat_t);
        
        free(dat_v);
        free(dat_t);
        
    }
  
    exit(EXIT_SUCCESS);
    
}

double test_func(double *p,double *k){
    return (cos(2*PI*k[0]*p[0])*cos(2*PI*k[1]*p[1])*cos(2*PI*k[2]*p[2]));
}


void writeSpline_test(int sType, string ofile, Triangulation *tri, double *dat_v, double *dat_t) {

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
    cout << "Writing to " << ofile << " nTtr="<<nTtr << std::endl;
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
    
    gfile << "</PointData>\n";
    
    gfile << "<CellData>\n";
    att.clear();
    att.push_back("type");  att.push_back("Float32");
    att.push_back("Name");  att.push_back("valueTtr");
    att.push_back("format");att.push_back("ascii");
    tri->TtrCentroidSplines(sType, dat_t);
    vtkXMLWriteDataArray(gfile,att,nTtr,dat_t);

    att[3] = "DxTtr";
    tri->TtrDerivativeSplines(sType,dat_t);
    vtkXMLWriteDataArray(gfile, att, nTtr, dat_t);
    
    
    
    gfile << "</CellData>\n";
    
    // FOOTER of VTK XML FILE
    vtkXMLFileClose(gfile);
    
    att.clear();
    free(idata);
    free(data);
}
