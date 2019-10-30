//
//  path.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 25.10.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "triangulate.hpp"
Path::Path(){
    return; 
}



Path::Path(int *p, int n,double tau_in,Triangulation *t){
    int i;
    vrt = std::vector<int>(p,p+n);
    if ( n == 0 )
        return; 
    edg.reserve(n-1);
    for ( i=0;i<vrt.size(); ++i)
        t->vrt[vrt[i]].nPth++;
    for (i=1;i<n; ++i)
        edg.push_back(t->isEdge(vrt[i],vrt[i-1]));
    tau=tau_in;
    return; 
}

int Path::len(){  return int(vrt.size()); }
