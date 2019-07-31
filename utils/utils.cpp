///
//  utils.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 16.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "utils.hpp"

template<class T> void sort4(T a[4])
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
template void sort4<int>   (int    a[4]);
template void sort4<double>(double a[4]);


template <class T>
int afind(T *v, int n, T val)
{
    for ( int i=0; i<n; ++i)
        if ( v[i] == val ) return(i);
    return n;
}
template int afind<string>(string *v, int n, string val);
template int afind<int>   (int    *v, int n, int val);
template int afind<double>(double *v, int n, double val);


template <class T>
int find(vector<T> vec,T val)
{
    for (int i=0;i<vec.size(); ++i)
        if (vec[i] == val)
            return i;
    return int(vec.size());
}
template int find<string>(vector<string>,string val);
template int find<double>(vector<double>,double val);
template int find<int>(vector<int>,int val);

uint64_t current_time() {
    using namespace std::chrono;
    return duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
   // return duration_cast<microseconds>(high_resolution_clock::now()).count();
}


bool file_exist(char *fc){
    string fs(fc);
    return file_exist(fs);
}
bool file_exist(string &fs){
    struct stat buffer;
    if ( stat (fs.c_str(), &buffer) != 0) {
        cout<<"Error: file \'" << fs << "\' not found \n";
        exit(EXIT_FAILURE);
    }
    return true;
}



void string_split(std::string& str, vector<string> &cont, char delim)
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}
