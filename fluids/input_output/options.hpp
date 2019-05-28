#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include "base_option.hpp"
#include <map>
#include <memory>
#include <string>
#include <vector>

class Options {
public:
    Options();
    Options(std::istream& is);

    void print_all();

    double gam(void) { return getDoubleOpt("GAMMA"); }
    double reynolds(void) { return getDoubleOpt("REYNOLDS"); }
    double prandtl(void) { return getDoubleOpt("PRANDTL"); }
    double mach(void) { return getDoubleOpt("MACH"); }
    double cfl(void) { return getDoubleOpt("CFL"); }
    double mix_param(void) { return getDoubleOpt("MIX_PARAM"); }
    double t_init(void) { return getDoubleOpt("T_INIT"); }
    double t_max(void) { return getDoubleOpt("T_MAX"); }
    double print_interval(void) { return getDoubleOpt("PRINT_INTERVAL"); }
    bool print_continue(void) { return getBoolOpt("PRINT_CONTINUE"); }

    std::string integrator_type(void)
    {
        return getStringOpt("INTEGRATOR_TYPE");
    }

    std::string solver_type(void) { return getStringOpt("SOLVER_TYPE"); }

    std::string input_base_path(void)
    {
        return getStringOpt("INPUT_BASE_PATH");
    }

    std::string input_meshfile_name(void)
    {
        return input_base_path() + getStringOpt("INPUT_MESHFILE_NAME");
    }

    std::string input_flow_configuration_file_name(void)
    {
        return input_base_path() + getStringOpt("INPUT_FLOW_CONFIGURATION");
    }

    std::string input_immersed_interface_file_name(void)
    {
        return input_base_path() + getStringOpt("INPUT_IMMERSED_INTERFACE");
    }

    std::string input_boundary_configuration_file_name(void)
    {
        return input_base_path() + getStringOpt("INPUT_BOUNDARY_CONFIGURATION");
    }

    std::string input_shock_file_name(void)
    {
        return input_base_path() + getStringOpt("INPUT_SHOCK_FILE");
    }

    std::string output_base_path(void)
    {
        return getStringOpt("OUTPUT_BASE_PATH");
    }
    std::string output_file_name(void)
    {
        return output_base_path() + getStringOpt("OUTPUT_FILE_NAME");
    }
    int output_counter(void) { return getIntOpt("OUTPUT_COUNTER"); }

    std::string input_type(void) { return getStringOpt("INPUT_TYPE"); }
    std::string output_type(void) { return getStringOpt("OUTPUT_TYPE"); }

    std::string flux(void) { return getStringOpt("FLUX"); }
    std::string convection(void) { return getStringOpt("CONVECTION"); }
    std::string mix_convection_main(void)
    {
        return getStringOpt("MIX_CONVECTION_MAIN");
    }
    std::string mix_convection_aux(void)
    {
        return getStringOpt("MIX_CONVECTION_AUX");
    }
    std::string dissipation(void) { return getStringOpt("DISSIPATION"); }

    int derivative_order(void) { return getIntOpt("DERIVATIVE_ORDER"); }

    int omp_threads(void) { return getIntOpt("OMP_THREADS"); }

    std::string luisa_detector() { return getStringOpt("LUISA_DETECTOR"); }

    double detector_sensitivity()
    {
        return getDoubleOpt("DETECTOR_SENSITIVITY");
    }

    bool should_filter() { return getBoolOpt("SHOULD_FILTER"); }
    int filter_order() { return getIntOpt("FILTER_ORDER"); }

private:
    std::map<std::string, std::unique_ptr<BaseOpt>> opt_map;

    double getDoubleOpt(const std::string& key);
    bool getBoolOpt(const std::string& key);
    std::string getStringOpt(const std::string& key);
    int getIntOpt(const std::string& key);

    std::map<std::string, std::unique_ptr<BaseOpt>> create_default_map(void);
    std::string get_input_base_path(std::istream& is);
    bool parse_file(std::istream& is);

    bool is_valid_opt(std::string& opt)
    {
        return (opt_map.find(opt) != opt_map.end());
    }

    bool safe_set_option(std::string& opt_name,
        std::vector<std::string>& opt_values, std::string& error_string);
    void got_invalid_option(std::string& opt_name, std::string& error_string);
};

#endif /* OPTIONS_HPP */
