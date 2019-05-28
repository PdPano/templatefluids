#include "tokenizer.hpp"
#include <algorithm>
#include <iostream>

inline void StringToUpperCase(std::string* str)
{
    std::transform((*str).begin(), (*str).end(), (*str).begin(),
        [](char c) { return (std::toupper(c)); });
}

bool TokenizeString(
    std::string& str, std::string& option_name, std::string& option_value)
{
    std::vector<std::string> wrapper = {option_value};
    bool ret = TokenizeString(str, option_name, wrapper);
    if (wrapper.size() != 1 and ret) {
        std::cerr << "Wrong number of elements returned! "
                  << "Expected 1, got " << wrapper.size() << std::endl;
        throw(-1);
    }
    option_value = wrapper[0];
    return ret;
}
bool TokenizeString(std::string& str, std::string& option_name,
    std::vector<std::string>& option_value)
{
    const std::string delimiters(" ()[]{}:,\t\n\v\f\r");
    // check for comments or empty string
    std::string::size_type pos, last_pos;
    pos = str.find_first_of('%');

    if ((str.length() == 0) || (pos == 0)) {
        // str is empty or a comment line, so no option here
        return false;
    }

    if (pos != std::string::npos) {
        // remove comment at end if necessary
        str.erase(pos);
    }

    // look for line composed on only delimiters (usually whitespace)
    pos = str.find_first_not_of(delimiters);

    if (pos == std::string::npos) {
        return false;
    }

    // find the equals sign and split string
    std::string name_part, value_part;
    pos = str.find('=');

    if (pos == std::string::npos) {
        std::cerr << "Error in TokenizeString(): "
                  << "line in the configuration file with no \"=\" sign."
                  << std::endl;
        std::cout << "Look for: " << str << std::endl;
        std::cout << "str.length() = " << str.length() << std::endl;
        throw(-1);
    }

    name_part = str.substr(0, pos);
    value_part = str.substr(pos + 1, std::string::npos);

    // the first_part should consist of one string with no interior delimiters
    last_pos = name_part.find_first_not_of(delimiters, 0);
    pos = name_part.find_first_of(delimiters, last_pos);

    if ((name_part.length() == 0) || (last_pos == std::string::npos)) {
        std::cerr << "Error in TokenizeString(): "
                  << "line in the configuration file with no name before the "
                     "\"=\" sign."
                  << std::endl;
        throw(-1);
    }

    if (pos == std::string::npos) {
        pos = name_part.length();
    }

    option_name = name_part.substr(last_pos, pos - last_pos);
    last_pos = name_part.find_first_not_of(delimiters, pos);

    if (last_pos != std::string::npos) {
        std::cerr << "Error in TokenizeString(): "
                  << "two or more options before an \"=\" sign in the "
                     "configuration file."
                  << std::endl;
        throw(-1);
    }

    StringToUpperCase(&option_name);

    // now fill the option value vector
    option_value.clear();
    last_pos = value_part.find_first_not_of(delimiters, 0);
    pos = value_part.find_first_of(delimiters, last_pos);

    while (std::string::npos != pos || std::string::npos != last_pos) {
        // add token to the vector<string>
        option_value.push_back(value_part.substr(last_pos, pos - last_pos));
        // skip delimiters
        last_pos = value_part.find_first_not_of(delimiters, pos);
        // find next "non-delimiter"
        pos = value_part.find_first_of(delimiters, last_pos);
    }

    if (option_value.empty()) {
        std::cerr << "Error in TokenizeString(): "
                  << "option " << option_name
                  << " in configuration file with no value assigned."
                  << std::endl;
        throw(-1);
    }

    // look for ';' DV delimiters attached to values
    std::vector<std::string>::iterator it;
    it = option_value.begin();

    while (it != option_value.end()) {
        if (*it == ";") {
            it++;
            continue;
        }

        pos = it->find(';');

        if (pos != std::string::npos) {
            std::string before_semi = it->substr(0, pos);
            std::string after_semi = it->substr(pos + 1, std::string::npos);

            if (before_semi.empty()) {
                *it = ";";
                it++;
                option_value.insert(it, after_semi);
            }
            else {
                *it = before_semi;
                it++;
                std::vector<std::string> to_insert;
                to_insert.emplace_back(";");

                if (!after_semi.empty()) {
                    to_insert.push_back(after_semi);
                }

                option_value.insert(it, to_insert.begin(), to_insert.end());
            }

            it = option_value.begin(); // go back to beginning; not efficient
            continue;
        }
        it++;
    }

    // remove any consecutive ";"
    it = option_value.begin();
    bool semi_at_prev = false;

    while (it != option_value.end()) {
        if (semi_at_prev) {
            if (*it == ";") {
                option_value.erase(it);
                it = option_value.begin();
                semi_at_prev = false;
                continue;
            }
        }

        semi_at_prev = *it == ";";

        it++;
    }

    return true;
}
