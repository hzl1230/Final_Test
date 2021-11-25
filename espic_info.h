#ifndef ESPIC_INFO_H
#define ESPIC_INFO_H

#include <string>

void espic_error(const std::string&);
void espic_warning(const std::string&);
std::string unknown_cmd_info(const std::string&, const std::string&);
std::string illegal_cmd_info(const std::string&, const std::string&);
std::string unknown_key_info(const std::string&, const std::string&, const std::string&);
std::string out_bound_info(const std::string&, const std::string& );
std::string illegal_arg_key_info(const std::string&, const std::string&, const std::string&, const std::string&);

#endif
