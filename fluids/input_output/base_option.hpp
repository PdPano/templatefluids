#ifndef BASE_OPTION_HPP
#define BASE_OPTION_HPP

#include <string>
struct BaseOpt {
    virtual void set(const std::string& in) = 0;
    virtual ~BaseOpt(){};
    virtual std::string print() const { return std::string("no"); };
};

#endif /* BASE_OPTION_HPP */
