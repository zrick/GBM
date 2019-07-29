///
//  utils.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 16.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "utils.hpp"

void sort4(int a[4])
{
    if (a[0]>a[1]) swap(a[0],a[1]);  // first two values in order
    if (a[1]>a[2]) {
        swap(a[1],a[2]);
        if (a[0]>a[1]) swap(a[0],a[1]);
    }                                // first three values in order
    if (a[2]>a[3]) {
        swap(a[2],a[3]);
        if (a[1]>a[2]) swap(a[1],a[2]);
        if (a[0]>a[1]) swap(a[0],a[1]);
    }                                // all (four) values in order
    
    return;
}


int ifind(int *v, int n, int val)
{
    for ( int i=0; i<n; ++i)
        if ( v[i] == val ) return(i);
    return n;
}

uint64_t current_time() {
    using namespace std::chrono;
    return duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
   // return duration_cast<microseconds>(high_resolution_clock::now()).count();
}

