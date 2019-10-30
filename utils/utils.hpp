//
//  utils.hpp
//  utils
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef utils_
#define utils_

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <stdio.h>
#include <sys/stat.h>
#include <utility>
#include <vector>

#include "constants.h"
#include "errors.hpp"

using namespace std;
void * GBMMalloc (int n);
void * GBMRealloc(void *p,int n);


void GBMError(string loc, string msg, int stat);
void GBMWarning(string loc, string msg, int stat);
void GBMLog(string msg);
void GBMOut(string msg);
void GBMASCWrite(string fname, string msg); 

void hashBoxIndices(double *p, double *dx, double *i);

template<class T>int afind(T *v, int n, T val);
template<class T>int find(vector<T> vec,T val);
template<class T>void sort4(T a[4]);
template<class T>int minpos(T *a, int n); 
uint64_t current_time();
bool file_exist(char *fc);
bool file_exist(string &fs);

void getCSVSize(string fname, int &ncol, int &nrow, char sep=GBM_DEFAULT_SEPARATOR);
void readCSV(string fname, vector<string> &header, vector<string> &labels,
             int &n, double *bare, char sep=GBM_DEFAULT_SEPARATOR);
void parseCSVLine(string line,  int iLine, double *v, string &lineLabel, int n=-1, char sep=GBM_DEFAULT_SEPARATOR);
void readCSVLine (string fname, int iLine, double *v, string &lineLabel, int n=-1, char sep=GBM_DEFAULT_SEPARATOR);

int string_split(std::string &str, vector<string> &cont, char delim=' ');

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ");

std::string gbmTime(); 



#endif
