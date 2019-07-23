//
//  main.cpp
//  vlinear
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include <iostream>
#include "linear.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    
    int nn=0,sing=0;
    double a[3],b[3],c[3];
    cout <<"Testing QR decomposition\n";
    
    sing=qrdcmp(nn);
    
    return 0;
}
