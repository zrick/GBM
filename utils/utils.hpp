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
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <stdio.h>
#include <sys/stat.h>
#include <utility>
#include <vector>

#include "constants.h"

using namespace std;

template<class T>int afind(T *v, int n, T val);
template<class T>int find(vector<T> vec,T val);
template<class T>void sort4(T a[4]);
template<class T>int minpos(T *a, int n); 
uint64_t current_time();
bool file_exist(char *fc);
bool file_exist(string &fs);

void string_split(std::string &str, vector<string> &cont, char delim=' ');

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ");
#endif
