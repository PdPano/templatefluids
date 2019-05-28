#ifndef DEFAULT_READER_HPP
#define DEFAULT_READER_HPP

#include "../options.hpp"

#include <fstream>
#include <iostream>
#include <vector>

#include "../../utils/boundary_point_def.hpp"
#include "../../utils/grid_constants_container_def.hpp"

struct GridComponentsContainer;

void default_reader(Options &opt, GridComponentsContainer &components,
                    GridConstantsContainer &constants);

void default_reader(Options &opt, GridComponentsContainer &components,
                    GridConstantsContainer &constants,
                    std::istream &&initial_conditions,
                    std::istream &&mesh_details, std::istream &&boundary_file,
                    std::istream &&immersed_interface,
                    std::istream &&shock_file);

GridConstantsContainer read_details_constants(std::istream &mesh_details);
std::pair<int, int> read_size_initial(std::istream &initial);
bool sizes_are_equal(std::pair<int, int> mesh, std::pair<int, int> initial);
std::vector<BoundaryPoint> read_boundary(std::istream &boundary_file);
void read_ghias_interface(GridComponentsContainer &components,
                          GridConstantsContainer &constants,
                          std::istream &immersed_interface_file);
void read_karagiozis_interface(GridComponentsContainer &components,
                               GridConstantsContainer &constants,
                               std::istream &immersed_interface_file);

void read_immersed_interface(GridComponentsContainer &components,
                             GridConstantsContainer &constants,
                             std::istream &immersed_interface);
void read_shock_points(GridComponentsContainer &components,
                       GridConstantsContainer &constants,
                       std::istream &shock_file, Options &opt);
#endif /* DEFAULT_READER_HPP */
