//
//  namelist.cpp
//  io
//
//  Created by Cedrick Ansorge on 14.08.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "io.hpp"

Namelist::Namelist(){
    return;
}

Namelist::Namelist(string &fname){
    read_namelist(fname);
    return; 
}

Namelist::Namelist(const char fname[]){
    string fstr(fname);
    read_namelist(fstr);
    return;
}

Namelist::Namelist(char *fname){
    string fstr(fname);
    read_namelist(fstr);
    return;
}

void Namelist::read_namelist(string s){
    int iLine=0;
    string line;
    
    if ( file_exist(s) ) {
        cout << "Reading Namelist from file " << s << "\n";
    }
    
    ifstream tfile(s);
    
    while ( getline(tfile,line) ) {
        iLine++;
        line=trim(line);
        // Skip comments
        if ( line[0] == '#' or line[0] == '!' or line[0] == '%' ) {
            continue;
        } else if ( line[0] == '[') {
            line=trim(line.erase(line.find_first_of(']')).erase(0,1));
            // Raise an error if the group already exists
            if ( find(grpNames,line) < grpNames.size() ) {
                cout << "ERROR: Group \'" << line << "\' exists alread\n";
                exit(EXIT_FAILURE);
            }
            grp.push_back(Group(line));
            grpNames.push_back(line); 
        } else if ( line.size() > 0 ) {
            std::vector<string> LineSplit;
            string_split(line,LineSplit,'=');
            if ( LineSplit.size() == 2 )
                grp[grp.size()-1].addAttribute(LineSplit[0],LineSplit[1]);
        }
    }
    
    /*
    for( i=0; i<grp.size(); ++i)
        for ( j=0; j<grp[i].att.size(); ++j)
            cout << grpNames[i] << ":" << grp[i].att[j].attName<< '=' << grp[i].att[j].attVal <<'\n';
    */

    return;
}

bool Namelist::hasVal(string group, string var){
    int iGrp;
    iGrp=find(grpNames,group);
    
    if ( iGrp >= grpNames.size() )
        return false;
    else {
        if ( grp[iGrp].hasAttribute(var) )
            return true;
        else
            return false;
    }
}

string Namelist::getValStr(string group,string var){
    int iGrp;
    Group *p_grp;
    string val_str;
    
    iGrp=find(grpNames,group);
    if ( iGrp >= grpNames.size() ) {
        cout << "ERROR: Group \'" << group << "\' not found in grpNames of len " << grpNames.size() << "\n";
        exit(EXIT_FAILURE);
    }
    p_grp=&grp[iGrp];
    p_grp->getAttribute(var, val_str);

    return val_str;
}

int Namelist::getVal_int(string group,string var) {
    return atoi(getValStr(group,var).c_str());
}

void Namelist::getList_dbl(string group, string var,double *list, const int nmax){
    int i;
    string s=getValStr(group,var);
    vector<string> s_split;
    string_split(s, s_split,',');
    if ( s_split.size() > nmax ){
        cout << "Warning: Found more values in \'" << var <<"\' of group \'" << group << "\' than needed\n";
        cout << "         Ignoring " << s_split.size()-nmax << "values.\n";
    }
    
    for ( i=0; i<nmax; ++i)
        list[i]=atof(s_split[i].c_str());
    
    return;
}


void Namelist::getList_int(string group, string var,int *list, const int nmax){
    cout << "not implemented \n";
    exit(EXIT_FAILURE);
    return; 
}

double Namelist::getVal_dbl(string group,string var) {
    return atof(getValStr(group,var).c_str());
}
                 
string Namelist::getVal_str(string group,string var) {
    return getValStr(group,var);
}

bool Namelist::getVal_bool(string group, string var) {
    string s=getVal_str(group,var);
    
    transform(s.begin(),s.end(),s.begin(),::tolower);

    if ( s.compare("t")==0 || s.compare("true")==0 || s.compare("0")==0 || s.compare(".true.")==0 )
        return true;
    if ( s.compare("f")==0 || s.compare("false")==0|| s.compare("1")==0 || s.compare(".false.")==0 )
        return false;
    else {
        cout << "ERROR: cannot guess bool from \'" << s << "\' for variable \'" << var <<"\' \n";
        exit(EXIT_FAILURE);
    }
}
