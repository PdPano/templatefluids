#include "stream_from_file.hpp"
#include <iostream>
std::ifstream stream_from_file(const std::string& filename)
{
    std::ifstream filestream;
    filestream.open(filename);
    if (!filestream.is_open()) {
        std::cerr << "Failed to open " << filename << std::endl;
        throw(-1);
    }
    return filestream;
}
