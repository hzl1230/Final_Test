#include <iostream>
#include <sstream>
#include <cstdlib>

#include "espic_info.h"

void espic_error(const std::string& msg)
{
  std::cerr << "--> Error: " << msg << ". <--" << std::endl;
  exit(-1);
}

void espic_warning(const std::string& msg)
{
  std::cerr << "--> Warning: " << msg << ". <--" << std::endl;
}

std::string unknown_cmd_info(const std::string& cmd, 
                             const std::string& file)
{
  std::ostringstream oss;
  oss << "Unknown command \"" << cmd << "\" in [" << file << "]";
  return oss.str();
}

std::string illegal_cmd_info(const std::string& cmd, 
                             const std::string& file)
{
  std::ostringstream oss;
  oss << "Illegal \"" << cmd << "\" command in [" << file << "]";
  return oss.str();
}

std::string out_bound_info(const std::string& var,
                           const std::string& file)
{
    std::ostringstream oss;
    oss << "Out the boundary of "<< var << " define in file [" << file << "]";
    return oss.str(); 
}


std::string unknown_key_info(const std::string& keyword, 
                             const std::string& cmd, 
                             const std::string& file)
{
  std::ostringstream oss;
  oss << "Unknown keyword \"" << keyword
      << "\" with command \"" << cmd << "\" in [" << file << "]";
  return oss.str();
}

std::string illegal_arg_key_info(const std::string& arg, const std::string& keyword,
                                 const std::string& cmd, const std::string& file)
{
  std::ostringstream oss;
  if (arg.empty()) {
    oss << "Illegal argument format in keyword \"" << keyword
      << "\" with command \"" << cmd << "\" in ["
      << file << "]";
  }
  else {
    oss << "Illegal argument \"" << arg
        << "\" of keyword \"" << keyword
        << "\" with command \"" << cmd << "\" in ["
        << file << "]";
  }

  return oss.str();
}