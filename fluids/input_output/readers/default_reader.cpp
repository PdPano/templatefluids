#include "default_reader.hpp"
#include "../../utils/body_discontinuity_def.hpp"
#include "../../utils/boundary_point_def.hpp"
#include "../../utils/generic_discontinuity_def.hpp"
#include "../../utils/ghias_ghost_point_def.hpp"
#include "../../utils/grid_components_container_def.hpp"
#include "../../utils/grid_constants_container_def.hpp"
#include "../../utils/point_def.hpp"
#include "../../utils/point_functions.hpp"
#include "../../utils/shock_discontinuity_handler.hpp"
#include "../stream_from_file.hpp"

void default_reader(Options& opt, GridComponentsContainer& components,
    GridConstantsContainer& constants)
{
    std::ifstream immersed_interface;
    std::ifstream shock_file;
    try {
        immersed_interface
            = stream_from_file(opt.input_immersed_interface_file_name());
    }
    catch (...) {
        std::cerr << "Immersed Interface File not found! Assuming empty!!!"
                  << std::endl;
    }

    if (opt.solver_type() == "SHOCK") {
        try {
            shock_file = stream_from_file(opt.input_shock_file_name());
        }
        catch (...) {
            std::cerr << "Shock file not found! Assuming empty!!!" << std::endl;
        }
    }
    default_reader(opt, components, constants,
        std::move(stream_from_file(opt.input_flow_configuration_file_name())),
        std::move(stream_from_file(opt.input_meshfile_name())),
        std::move(stream_from_file(
            opt.input_boundary_configuration_file_name())),
        std::move(immersed_interface), std::move(shock_file));
}

void default_reader(Options& opt, GridComponentsContainer& components,
    GridConstantsContainer& constants, std::istream&& initial_conditions,
    std::istream&& mesh_details, std::istream&& boundary_file,
    std::istream&& immersed_interface, std::istream&& shock_file)
{
    constants = read_details_constants(mesh_details);

    std::pair<int, int> mesh_size, initial_size;
    mesh_size = std::make_pair(constants.nPointsI, constants.nPointsJ);

    initial_size = read_size_initial(initial_conditions);

    if (!sizes_are_equal(mesh_size, initial_size)) {
        std::cerr << "Input mesh file: " << opt.input_meshfile_name()
                  << std::endl
                  << "Input initial conditions: "
                  << opt.input_flow_configuration_file_name() << std::endl;
    }

    components.grid_c = std::vector<Point>(constants.nPointsTotal);
    components.flags_c = std::vector<int>(constants.nPointsTotal);

    for (int i = 0; i < constants.nPointsTotal; i++) {
        mesh_details >> components.flags_c[i];
        initial_conditions >> components.grid_c[i];
    }

    components.boundary_c = read_boundary(boundary_file);

    read_immersed_interface(components, constants, immersed_interface);

    read_shock_points(components, constants, shock_file, opt);
}

GridConstantsContainer read_details_constants(std::istream& mesh_details)
{
    GridConstantsContainer constants{};
    int nPointsI, nPointsJ, nPointsTotal;
    double dx, dy, xmin, ymin;

    mesh_details >> nPointsI >> nPointsJ;
    nPointsTotal = nPointsI * nPointsJ;

    mesh_details >> dx >> dy >> xmin >> ymin;

    constants.nPointsI = nPointsI;
    constants.nPointsJ = nPointsJ;
    constants.nPointsTotal = nPointsTotal;
    constants.dx = dx;
    constants.dy = dy;
    constants.xmin = xmin;
    constants.ymin = ymin;

    return constants;
}

std::pair<int, int> read_size_initial(std::istream& initial)
{
    int nPointsI, nPointsJ;
    initial >> nPointsI >> nPointsJ;
    return std::make_pair(nPointsI, nPointsJ);
}

bool sizes_are_equal(std::pair<int, int> mesh, std::pair<int, int> initial)
{
    if (mesh != initial) {
        std::cerr << "Sizes don't match!" << std::endl
                  << "Read from mesh: " << mesh.first << " " << mesh.second
                  << std::endl
                  << "Read from initial condition: " << initial.first << " "
                  << initial.second << std::endl;
        return false;
    }
    return true;
}

std::vector<BoundaryPoint> read_boundary(std::istream& boundary_file)
{
    int num_points;
    boundary_file >> num_points;
    std::vector<BoundaryPoint> boundary(num_points);
    for (int i = 0; i < num_points; i++) {
        boundary_file >> boundary[i];
    }
    return boundary;
}

void read_immersed_interface(GridComponentsContainer& components,
    GridConstantsContainer& constants, std::istream& immersed_interface)
{
    std::string interface_type;
    immersed_interface >> interface_type;
    if (interface_type == "GHIAS") {
        read_ghias_interface(components, constants, immersed_interface);
    }
    else if (interface_type == "KARAGIOZIS") {
        read_karagiozis_interface(components, constants, immersed_interface);
    }
}

void read_ghias_interface(GridComponentsContainer& components,
    GridConstantsContainer& constants, std::istream& immersed_interface_file)
{
    auto get_ind = [&](int i, int j) { return i * constants.nPointsJ + j; };

    int num_points;
    immersed_interface_file >> num_points;
    components.ghias_points_c.reserve(num_points);

    for (int gp = 0; gp < num_points; gp++) {
        GhiasGhostPoint ghost{};
        int i, j;
        double x, y, nx, ny;
        char type;

        immersed_interface_file >> i >> j;
        ghost.ind = get_ind(i, j);
        for (int k = 0; k < 4; k++) {
            immersed_interface_file >> type >> i >> j >> x >> y >> nx >> ny;
            ghost.is_fluid_point[k] = (type == 'f'); // NOLINT
            ghost.neighbors_inds[k] = get_ind(i, j); // NOLINT
            ghost.neighbors_x[k] = x;                // NOLINT
            ghost.neighbors_y[k] = y;                // NOLINT
            ghost.neighbors_nx[k] = nx;              // NOLINT
            ghost.neighbors_ny[k] = ny;              // NOLINT
        }
        immersed_interface_file >> x >> y;
        ghost.image_coordinate[0] = x;
        ghost.image_coordinate[1] = y;
        components.ghias_points_c.push_back(ghost);
    }
}

void read_karagiozis_interface(GridComponentsContainer& components,
    GridConstantsContainer& constants, std::istream& immersed_interface_file)
{
    auto get_ind = [&](int i, int j) { return i * constants.nPointsJ + j; };

    int num_points;
    immersed_interface_file >> num_points;
    components.karagiozis_points_c.reserve(num_points);

    for (int bp = 0; bp < num_points; bp++) {
        BodyDiscontinuity bd{};
        int i, j;

        immersed_interface_file >> bd.type;
        immersed_interface_file >> i >> j;
        bd.ind = get_ind(i, j);
        immersed_interface_file >> bd.frac;
        immersed_interface_file >> bd.theta;
        immersed_interface_file >> bd.left() >> bd.right();
        components.karagiozis_points_c.push_back(bd);
    }
}

void read_shock_points(GridComponentsContainer& components,
    GridConstantsContainer& constants, std::istream& shock_file, Options& opt)
{
    std::string data_type;
    shock_file >> data_type;
    if (data_type == "SHOCK") {
        ShockHandler sh(opt);
        auto get_ind = [&](int i, int j) { return i * constants.nPointsJ + j; };
        int num_points;
        shock_file >> num_points;
        components.shock_points_c.reserve(num_points);
        char type;
        int ind;
        double frac;
        Point left_p, right_p;

        for (int sp = 0; sp < num_points; sp++) {
            int i, j;

            shock_file >> type;
            shock_file >> i >> j;
            ind = get_ind(i, j);
            shock_file >> frac;
            shock_file >> left_p >> right_p;
            auto sd = sh.create_shock(left_p, right_p, type);
            sd.frac = frac;
            sd.ind = ind;
            components.shock_points_c.push_back(sd);
        }
    }
}
