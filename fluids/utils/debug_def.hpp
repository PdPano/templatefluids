#ifndef DEBUG_DEF_HPP
#define DEBUG_DEF_HPP

#ifdef DEBUG
#include <iostream>
#define DUMP(x) std::cout << #x << "=" << (x) << std::endl
#else
#define DUMP(x)
#endif

#endif /* DEBUG_DEF_HPP */
