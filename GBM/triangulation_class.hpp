//
//  triangulation_class.hpp
//  GBM
//
//  Created by Cedrick Ansorge on 11.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef triangulation_class_hpp

class Triangulation{
    
public:
    // Constructors
    Triangulation();
    Triangulation(char *fname);
    
    // Members
    int status;
    
private:
    int parse_triangle_file(char *fname);
    
};

#define triangulation_class_hpp

#include <stdio.h>

#endif /* triangulation_class_hpp */
