//
//  io.cpp
//  io
//
//  Created by Cedrick Ansorge on 30.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "io.hpp"

Attribute::Attribute(string typ, string val) {
    typ=string(typ);
    val=string(val);
    return;
}

void vtkXMLFileOpen(string &name,int np,int nc,ofstream &f) {
    f=ofstream(name);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<UnstructuredGrid>\n";
    f << "<Piece NumberOfPoints=\"" << np <<"\" NumberOfCells=\"" << nc << "\">\n";
    return;
}

void vtkXMLFileClose(ofstream &f){
    f << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    return;
}

template<typename T>
void vtkXMLWriteDataArray(ofstream &gfile, vector<string> att, int nval, T *val){
    int i;
    const int natt = int(att.size()/2);
    
    gfile << "<DataArray ";
    for (i=0;i<natt;++i)
        gfile << att[2*i] << "=\"" << att[2*i+1] << "\" ";
    gfile << ">\n";
    for (i=0;i<nval;++i)
        gfile << val[i] << " "; 
    gfile << "\n</DataArray>\n";
    
   
    return;
}
template void vtkXMLWriteDataArray<double>(ofstream &gfile, vector<string> att, int nval, double *val);
template void vtkXMLWriteDataArray<int>(ofstream &gfile, vector<string> att, int nval, int *val);

