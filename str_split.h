#ifndef STR_SPLIT_H
#define STR_SPLIT_H

#include <string>
#include <sstream>
#include <vector>
#include <cstring>

const int MaxLine = 1024*1024;

using std::string;
using std::stringstream;
using std::vector;

void split(const string& , char , vector<string>& );
vector<string> split(const string& , char );

void split(const string& , char* , vector<string>& );
vector<string> split(const string& , char* );

#endif
