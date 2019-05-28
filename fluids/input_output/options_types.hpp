#ifndef OPTIONS_TYPES_HPP
#define OPTIONS_TYPES_HPP

#include "base_option.hpp"
#include <algorithm>
#include <sstream>
#include <string>

struct IntOpt : public BaseOpt {

    IntOpt(const std::string& in)
        : val(parse(in)){};

    void set(const std::string& in) final { val = parse(in); };
    int value(void) const { return val; }
    std::string print() const { return std::to_string(val); }

private:
    int val;
    int parse(const std::string& in) { return std::stoi(in); }
};

struct BoolOpt : public BaseOpt {

    BoolOpt(const std::string& in)
        : val(parse(in)){};

    void set(const std::string& in) final { val = parse(in); };
    bool value(void) const { return val; }
    std::string print() const { return (val) ? "TRUE" : "FALSE"; }

private:
    bool val;
    bool parse(const std::string& in)
    {
        bool ret;
        std::string local;
        std::transform(
            in.begin(), in.end(), std::back_inserter(local), ::tolower);
        std::istringstream(local) >> std::boolalpha >> ret;
        return ret;
    }
};

struct DoubleOpt : public BaseOpt {

    DoubleOpt(const std::string& in)
        : val(parse(in)){};

    void set(const std::string& in) final { val = parse(in); };
    double value(void) const { return val; }
    std::string print() const { return std::to_string(val); }

private:
    double val;
    double parse(const std::string& in) { return std::stod(in); }
};

struct StringOpt : public BaseOpt {

    StringOpt(std::string in)
        : val(in){};

    void set(const std::string& in) final { val = in; };
    std::string value(void) const { return val; }
    std::string print() const { return val; }

private:
    std::string val;
};

#endif /*OPTIONS_TYPES_HPP*/
