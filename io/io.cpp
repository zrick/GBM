//
//  io.cpp
//  io
//
//  Created by Cedrick Ansorge on 30.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "io.hpp"

Group::Group(string &s){
    name=string(s);
    
    return;
}

void Group::addAttribute(string n, string v) {
    int i;
    
    for (i=0; i<att.size(); ++i )
        if ( att[i].attName.compare(n) == 0 ){
            cout << "ERROR: Value for variable \'" << n <<"\' in group \'" << name << "\' exists already\n";
            exit(EXIT_FAILURE);
        }
    
    att.push_back(Attribute(n,v));
    
    return;
}

bool Group::hasAttribute(string &n) {
    for (int i=0; i<att.size(); ++i )
        if ( att[i].attName.compare(n) == 0 ){
            return true;
        }
    return false;
}

void Group::getAttribute(string &n, string &v) {
    int i;
    
    for (i=0; i<att.size(); ++i )
        if ( att[i].attName.compare(n) == 0 ){
            v = string(att[i].attVal);
            return;
        }
    
    if ( i == att.size() )
    {
        cout << "Error: Variable \'" << n << "\' not found in group \'" << name << "\' with " << att.size() << " variables\n";
        exit(EXIT_FAILURE);
    }
    
    return;
}

Attribute::Attribute(string &typ, string &val) {
    attName=string(typ);
    attVal =string(val);
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
    f.close();
    return;
}

template<typename T>
void vtkXMLWriteDataArray(ofstream &gfile, vector<string> att, int nval, T *val){
    int i;
    const int natt = int(att.size()/2);
    char *buf;
    stringstream st;
    int nchar=0;
    string format;
    
    buf=(char *)malloc((nval*10+100)*sizeof(char)); 
    
    if (typeid(T) == typeid(int) )
        format = string("%d ");
    else if ( typeid(T) == typeid(double) )
        format = string("%9.5f ");
    
    nchar += sprintf(&buf[nchar],"<DataArray ");
    for (i=0;i<natt;++i)
        nchar += sprintf(&buf[nchar],"%s=\"%s\" ", att[2*i].c_str(), att[2*i+1].c_str());
    nchar += sprintf(&buf[nchar],">\n");
    for (i=0;i<nval;i++)
        nchar += sprintf(&buf[nchar],format.c_str(),val[i]);
    nchar += sprintf(&buf[nchar],"\n</DataArray>\n");
    gfile.write(buf,nchar);
   
    return;
}
template void vtkXMLWriteDataArray<double>(ofstream &gfile, vector<string> att, int nval, double *val);
template void vtkXMLWriteDataArray<int>(ofstream &gfile, vector<string> att, int nval, int *val);

