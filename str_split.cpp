#include "str_split.h"

void split(const string &s, char delim, vector<string> &elems) 
{
  stringstream ss;
  ss.str(s);
  string item;
  while (std::getline(ss, item, delim)) {
    if (!item.empty())
      elems.push_back(item);
  }
}

vector<string> split(const string &s, char delim) 
{
  vector<string> elems;
  split(s, delim, elems);
  return elems;
}

void split(const string &s, char* delim, vector<string> &elems) 
{
  char cstr[MaxLine];
  memcpy(cstr, s.data(), s.length());
  cstr[s.length()] = '\0';
//  for (uint32_t i = 0; i <= s.length(); ++i) cstr[i] = s[i];

  char *token = std::strtok(cstr, delim);
  while (token) {
    elems.push_back(string(token));
    token = std::strtok(NULL, delim);
  }
}

vector<string> split(const string &s, char* delim) 
{
  vector<string> elems;
  split(s, delim, elems);
  return elems;
}
