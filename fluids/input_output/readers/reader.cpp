#include "reader.hpp"

Reader::Reader(Options& opt)
    : local_container(GridConstantsContainer())
{
    std::string input_type = opt.input_type();
    if (input_type == "DEFAULT") {
        default_reader(opt, local_components, local_container);
    }
    else {
        std::cerr << "Input type '" << input_type << "' is not supported"
                  << std::endl;
        throw(-1);
    }
}

Reader::Reader(Options& opt, std::istream&& initial_conditions,
    std::istream&& mesh_details, std::istream&& boundary_file,
    std::istream&& immersed_interface, std::istream&& shock_file)
    : local_container(GridConstantsContainer())
{
    std::string input_type = opt.input_type();
    if (input_type == "DEFAULT") {
        default_reader(opt, local_components, local_container,
            std::move(initial_conditions), std::move(mesh_details),
            std::move(boundary_file), std::move(immersed_interface),
            std::move(shock_file));
    }
    else {
        std::cerr << "Input type '" << input_type << "' is not supported"
                  << std::endl;
        throw(-1);
    }
}
