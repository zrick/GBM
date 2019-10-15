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
    
    if ( file_exist(s) )
        GBMLog("Reading Namelist from file " + s);
    
    ifstream tfile(s);
    
    while ( getline(tfile,line) ) {
        iLine++;
        line=trim(line);
        // Skip comments
        if ( line[0] == '#' or line[0] == '!' or line[0] == '%' ) {
            continue;
        } else if ( line[0] == '[') {
            line=trim(line.erase(line.find_first_of(']')).erase(0,1));

            if ( find(grpNames,line) < grpNames.size() )
                GBMError("Namelist::read_namelist", "Group \'" + line + "\' exists already",GBMERR_NAMELIST);

            grp.push_back(Group(line));
            grpNames.push_back(line); 
        } else if ( line.size() > 0 ) {
            std::vector<string> LineSplit;
            string_split(line,LineSplit,'=');
            if ( LineSplit.size() == 2 )
                grp[grp.size()-1].addAttribute(LineSplit[0],LineSplit[1]);
        }
    }
    
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
    if ( iGrp >= grpNames.size() )
        GBMError("Namelist::getValStr", "Group \'" + group + "\' not found in grpNames of len " + to_string( grpNames.size()), GBMERR_NAMELIST);

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
    if ( s_split.size() > nmax ) {
        GBMWarning("Namelist::getList_dbl", "Found more values in \'" + var +"\' of group \'" + group + "\' than needed",GBMERR_NAMELIST);
        GBMWarning("Namelist::getList_dbl", "Ignoring " +to_string(s_split.size()-nmax) + "values.",GBMERR_NAMELIST);
    }
    
    for ( i=0; i<nmax; ++i)
        list[i]=atof(s_split[i].c_str());
    
    return;
}


void Namelist::getList_int(string group, string var,int *list, const int nmax){
    GBMError("Namelist:getList_int","Not Implemented",GBMERR_UNDEVELOPED);
    return; 
}

double Namelist::getVal_dbl(string group,string var) {
    if ( ! hasVal(group,var) )
        GBMError("Namelist:getVal_dbl", string("No default for [" + group +"]:" + var), GBMERR_NAMELIST);
    return atof(getValStr(group,var).c_str());
}
                 
string Namelist::getVal_str(string group,string var) {
    if ( !hasVal(group,var) )
        GBMError("Namelist:getVal_dbl", string("No default for [" + group +"]:" + var), GBMERR_NAMELIST);
    return getValStr(group,var);
}

bool Namelist::getVal_bool(string group, string var) {
    string s=getVal_str(group,var);

    if ( !hasVal(group,var) )
        GBMError("Namelist:getVal_dbl", string("No default for [" + group +"]:" + var), GBMERR_NAMELIST);
    
    transform(s.begin(),s.end(),s.begin(),::tolower);

    if ( s.compare("t")==0 || s.compare("true")==0 || s.compare("0")==0 || s.compare(".true.")==0 )
        return true;
    if ( s.compare("f")==0 || s.compare("false")==0|| s.compare("1")==0 || s.compare(".false.")==0 )
        return false;
    else {
        GBMError("Namelist::getVal_bool", "cannot guess bool from \'" + s + "\' for variable \'" + var  +"\'",GBMERR_NAMELIST);
        exit(EXIT_FAILURE);
    }
}
