//
//  triangulation_class.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "triangulation_class.hpp"
#include "main.hpp"


Triangulation::Triangulation(){
}

Triangulation::Triangulation(char *fname){
    status=parse_triangle_file(fname);
}

int Triangulation::parse_triangle_file(char *fname){
    
    using namespace std; 
    // LOCAL DECLARATIONS
    double a,b,c;
    int iLine=0;
    struct stat buffer;
    string line;
    ifstream tfile(fname);
        
    // check if file exists
    if ( stat (fname, &buffer) != 0) {
        cout<<"Error: file \'" << &fname[0] << "\' not found \n";
    } else {
        cout << "Reading from file \'" << &fname[0] << "\'\n";
    }
        
    for ( string line; getline(tfile,line); ) {
        istringstream lss(line);
        lss >> a >> b >> c;
        cout << "LINE "<< iLine << ":  " << a << "; " << b << "; "<< c << '\n';
        iLine++;
    }
    return(0);
}
