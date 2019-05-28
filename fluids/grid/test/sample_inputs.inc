std::string  grid_info_sample= (
    "2 2\n"
    "0.01 0.02\n"
    "-1.1 -1.2\n"
    "0 1\n"
    "0 1\n"
    );

std::string initial_conditions_sample = (
    "2 2\n"
    "1.0 0.0 0.0 40.0\n"
    "2.0 0.0 0.0 40.0\n"
    "2.0 0.0 0.0 40.0\n"
    "1.0 0.0 0.0 40.0\n"
    );

std::string boundary_sample = (
    "1\n"
    "1\n"
    "1 0\n"
    "2 0\n"
    "1.0 2.0 3.0 4.0 5.0\n"
    "1 2 3\n"
    );

std::string immersed_interface_all_fluid = (
    "GHIAS\n"
    "1\n"
    "1 0\n"
    "f 0 0 0.0 0.0 0.0 0.0\n"
    "f 0 1 1.0 0.0 0.0 0.0\n"
    "f 1 0 0.0 1.0 0.0 0.0\n"
    "f 1 1 1.0 1.0 0.0 0.0\n"
    "0.5 0.5\n"
    );
std::string immersed_interface_one_border = (
    "GHIAS\n"
    "1\n"
    "1 0\n"
    "f 0 0 0.0 0.0 0.0 0.0\n"
    "f 0 1 1.0 0.0 0.0 0.0\n"
    "b 1 0 0.25 0.75 0.707106781186 -0.707106781186\n"
    "f 1 1 1.0 1.0 0.0 0.0\n"
    "0.5 0.5\n"
    );
