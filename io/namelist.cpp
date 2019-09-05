//
//  namelist.cpp
//  io
//
//  Created by Cedrick Ansorge on 14.08.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "io.hpp"

Group :: Group(string &s){
    name=string(s);
    return; 
}

Namelist::Namelist(){
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
            cout << "Starting Group with line ";
            line=trim(line.erase(line.find_first_of(']')).erase(0,1));
            // Raise an error if the group already exists
            
            grp.push_back(Group(line));
            grpNames.push_back(line); 
        } else {
            cout << "Parsing Line             ";
        }
        cout << "Line " << iLine << ":" << line << "\n";
    }
    
    return;
}
