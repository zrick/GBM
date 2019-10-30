///
//  utils.cpp
//  GBM
//
//  Created by Cedrick Ansorge on 16.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#include "utils.hpp"

void * GBMMalloc(int n){
    void *p;
    if  ( ( p  = (void *) malloc(n*sizeof(p)) ) == nullptr )
        GBMError("libutils:gbmMalloc", "Allocation of " + to_string(n*sizeof(p)) +" bytes failed", GBMERR_MALLOC);
    else
        GBMLog("libutils:gbmMalloc: Allocated memory chunk of " + to_string(n*sizeof(p)) + " Bytes");
    return p;
}

void * GBMRealloc(void *p,int n){
    if  ( ( p = (void *) realloc(p,n*sizeof(p)) ) == nullptr )
        GBMError("libutils:gbmRealloc", "Allocation of " + to_string(n*sizeof(p)) +" bytes failed", GBMERR_MALLOC);
    else
        GBMLog("libutils:gbmRealloc: Re-allocated Memory Chunk of " + to_string(n*sizeof(p)) + " Bytes");
    return p;
}

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

int string_split(std::string& str, vector<string> &cont, char delim)
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
    return int(cont.size()); 
}

template<typename T>
int minpos(T *arr, int n) {
    int i,pos=0;
    T   v=0;
    for (i=0;i<n;++i)
        if ( arr[i] < v ) {
            pos=i;
            v=arr[i];
        }
    return pos;
}
template int minpos<double>(double *arr, int n);
template int minpos<int>(int *arr, int n);

void getCSVSize(string fname, int &ncol, int &nrow,char sep){
    vector<string> dum;
    string line;

    if ( ! file_exist(fname) )
        GBMError("libutils:getCSVSize", "File " + fname + " not found", GBMERR_GENERAL );
    
    ifstream f(fname);
    if ( getline(f,line)) { ncol=string_split(line,dum,sep); }

    nrow=0;
    while(getline(f,line)) { ++nrow;}
    return;
}

void readCSV(string fname, vector<string> &header, vector<string> &labels,
             int &n, double *bare, char sep){
    string l, lab;
    int iL=0;
    ifstream f;
    int iPtr=0, nHead=0;
    
    if ( ! file_exist(fname) )
        GBMError("libutils:readCSV", "File " + fname + " not found", GBMERR_GENERAL );
    
    f=ifstream(fname);
    
    if ( getline(f,l) ) { nHead=string_split(l,header,sep); ++iL; }
    
    iPtr=0;
    while ( getline(f,l) ) {
        parseCSVLine(l,iL,&bare[iPtr],lab,nHead);
        ++iL; iPtr+=nHead;
        labels.push_back(lab);
    }
    n=iL-1;
    
    return;
}

void readCSVLine(string fname, int iLine, double *v, string &lineLabel, int n, char sep) {
    vector<string> split;
    int l_loc=0;
    ifstream f(fname);
    string line;
    
    while ( getline(f,line) && l_loc < iLine ) ++l_loc;
    if  ( l_loc != iLine )
        GBMError("libutils:readCSVLine",
                 "Problem finding Line " + to_string(iLine) + " in file " + fname,
                 GBMERR_INDEX);
    parseCSVLine(line,iLine,v,lineLabel);
    return;
}

void parseCSVLine(string line, int iLine, double *v, string &lineLabel, int n, char sep) {
    vector<string> split;
    int nval;
    
    if ( (nval=string_split(line, split, sep)) != n && n>0 )
        GBMWarning("libutils:readCSV",
                   "Number of elements in line (" + to_string(iLine) +") does not match header size of " + to_string(n), GBMERR_INDEX);

    lineLabel=split[0];
    for (int i=0; i<nval; ++i)
        try { v[i] = stod(split[i].c_str()); }
        catch ( const std::invalid_argument ) { v[i] = i==0? (iLine-1) : GBMMissval; }
        catch ( const std::out_of_range )     { v[i] = GBMMissval;    }
    
    return;
}

std::string& ltrim(std::string& str, const std::string& chars){
    str.erase(0, str.find_first_not_of(chars));
    return str;
}

std::string& rtrim(std::string& str, const std::string& chars){
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}

std::string& trim(std::string& str, const std::string& chars){
    return ltrim(rtrim(str, chars), chars);
}

std::string gbmTime(){
    time_t _tm =time(NULL );
    stringstream s;
    struct tm * curtime = localtime ( &_tm );
    s<<asctime(curtime);
    return s.str();
}

void GBMError(string loc, string msg, int stat){
    GBMWarning(loc,msg,stat);
    exit(stat);
}

void GBMWarning(string loc, string msg, int stat){
    string msg_loc = ">>> GBM Warning " +to_string(stat) + " in "+ loc + ": " + msg;
    cout << msg_loc << std::endl;        // Warnings and errors are written to all log files and stdout 
    GBMASCWrite(GBM_FILE_ERR,msg_loc);
    GBMASCWrite(GBM_FILE_LOG,msg_loc);
    GBMASCWrite(GBM_FILE_OUT,msg_loc);
    return;
}

void GBMLog(string msg){ GBMASCWrite(GBM_FILE_LOG, msg); }
void GBMOut(string msg){ GBMASCWrite(GBM_FILE_OUT, msg); }

void GBMASCWrite(string fname, string msg) {
    ofstream os(fname,ios_base::app);
    os << msg << std::endl;
    os.close();
    return;
}
