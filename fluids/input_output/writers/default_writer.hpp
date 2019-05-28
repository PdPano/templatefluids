#ifndef DEFAULT_WRITER_HPP
#define DEFAULT_WRITER_HPP

#include <string>
class CartesianGrid;
class KaragiozisGrid;
struct PointFunctions;

void default_writer(
    CartesianGrid& grid, std::string& file_name, PointFunctions& pf);

void default_writer(
    KaragiozisGrid& grid, std::ofstream& output, PointFunctions& pf);

#endif /* DEFAULT_WRITER_HPP */
