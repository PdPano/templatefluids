#ifndef TOKENIZER_HPP
#define TOKENIZER_HPP

#include <string>
#include <vector>
bool TokenizeString(
    std::string& str, std::string& option_name, std::string& option_value);

bool TokenizeString(std::string& str, std::string& option_name,
    std::vector<std::string>& option_value);

#endif
