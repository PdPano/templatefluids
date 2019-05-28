#include "options.hpp"
#include "options_types.hpp"
#include "tokenizer.hpp"
#include <iostream>

Options::Options() : opt_map(Options::create_default_map()) {}

Options::Options(std::istream &is) : Options() {
  if (parse_file(is)) {
    std::cerr << "Options constructor failed";
    throw(-1);
  }
}

void Options::print_all() {
  for (auto &opt : opt_map) {
    std::cout << opt.first << " = " << opt.second->print() << std::endl;
  }
}

std::map<std::string, std::unique_ptr<BaseOpt>> Options::create_default_map() {

  std::map<std::string, std::unique_ptr<BaseOpt>> def_map;
  def_map["GAMMA"] = std::make_unique<DoubleOpt>("1.4");
  def_map["REYNOLDS"] = std::make_unique<DoubleOpt>("10.0");
  def_map["PRANDTL"] = std::make_unique<DoubleOpt>("1.0");
  def_map["MACH"] = std::make_unique<DoubleOpt>("0.1");
  def_map["CFL"] = std::make_unique<DoubleOpt>("0.9");
  def_map["MIX_PARAM"] = std::make_unique<DoubleOpt>("0.1");

  def_map["T_MAX"] = std::make_unique<DoubleOpt>("1.0");
  def_map["T_INIT"] = std::make_unique<DoubleOpt>("0.0");
  def_map["PRINT_INTERVAL"] = std::make_unique<DoubleOpt>("0.1");
  def_map["PRINT_CONTINUE"] = std::make_unique<BoolOpt>("FALSE");

  def_map["INTEGRATOR_TYPE"] = std::make_unique<StringOpt>("EULER");
  def_map["SOLVER_TYPE"] = std::make_unique<StringOpt>("SIMPLE");

  def_map["INPUT_BASE_PATH"] = std::make_unique<StringOpt>("./input/");
  def_map["INPUT_MESHFILE_NAME"] = std::make_unique<StringOpt>("gridInfo.dat");
  def_map["INPUT_IMMERSED_INTERFACE"] =
      std::make_unique<StringOpt>("immersedInterface.dat");
  def_map["INPUT_FLOW_CONFIGURATION"] =
      std::make_unique<StringOpt>("initialCondition.dat");
  def_map["INPUT_BOUNDARY_CONFIGURATION"] =
      std::make_unique<StringOpt>("boundary.dat");
  def_map["INPUT_SHOCK_FILE"] = std::make_unique<StringOpt>("shock.dat");
  def_map["OUTPUT_FILE_NAME"] = std::make_unique<StringOpt>("result");
  def_map["OUTPUT_BASE_PATH"] = std::make_unique<StringOpt>("./output/");
  def_map["OUTPUT_COUNTER"] = std::make_unique<IntOpt>("0");
  def_map["INPUT_TYPE"] = std::make_unique<StringOpt>("DEFAULT");
  def_map["OUTPUT_TYPE"] = std::make_unique<StringOpt>("VTK");

  def_map["FLUX"] = std::make_unique<StringOpt>("SIMPLE");
  def_map["CONVECTION"] = std::make_unique<StringOpt>("SIMPLE");
  def_map["MIX_CONVECTION_MAIN"] = std::make_unique<StringOpt>("SIMPLE");
  def_map["MIX_CONVECTION_AUX"] =
      std::make_unique<StringOpt>("SPLIT_CONVECTION");
  def_map["DISSIPATION"] = std::make_unique<StringOpt>("SIMPLE");

  def_map["DERIVATIVE_ORDER"] = std::make_unique<IntOpt>("2");

  def_map["OMP_THREADS"] = std::make_unique<IntOpt>("0");

  def_map["LUISA_DETECTOR"] = std::make_unique<StringOpt>("TYPE_23");
  def_map["DETECTOR_SENSITIVITY"] = std::make_unique<DoubleOpt>("1.0");

  def_map["SHOULD_FILTER"] = std::make_unique<BoolOpt>("FALSE");
  def_map["FILTER_ORDER"] = std::make_unique<IntOpt>("3");

  return def_map;
}

double Options::getDoubleOpt(const std::string &key) {
  return dynamic_cast<DoubleOpt *>(opt_map[key].get())->value();
}

int Options::getIntOpt(const std::string &key) {
  return dynamic_cast<IntOpt *>(opt_map[key].get())->value();
}

bool Options::getBoolOpt(const std::string &key) {
  return dynamic_cast<BoolOpt *>(opt_map[key].get())->value();
}

std::string Options::getStringOpt(const std::string &key) {
  return dynamic_cast<StringOpt *>(opt_map[key].get())->value();
}

bool Options::parse_file(std::istream &is) {
  std::string base_path = get_input_base_path(is);
  if (base_path.length() != 0u) {
    opt_map["INPUT_BASE_PATH"]->set(base_path);
  }
  std::string line, opt_name;
  std::vector<std::string> opt_values;
  std::string error_string;
  int error_count = 0;
  while (getline(is, line)) {
    if (TokenizeString(line, opt_name, opt_values)) {
      if (opt_name[0] == '#') { // Comment line
        continue;
      }
      if (is_valid_opt(opt_name)) {
        if (!safe_set_option(opt_name, opt_values, error_string)) {
          error_count++;
        }
      } else {
        got_invalid_option(opt_name, error_string);
        error_count++;
      }
    }
  }
  if (!error_string.empty()) {
    std::cerr << error_string << std::endl;
    return true;
  }
  return false;
}

bool Options::safe_set_option(std::string &opt_name,
                              std::vector<std::string> &opt_values,
                              std::string &error_string) {
  try {
    opt_map[opt_name]->set(opt_values[0]);
    return true;
  } catch (...) {
    std::string error_msg;
    error_msg.append("Something wrong during assingment ");
    error_msg.append(opt_name);
    error_msg.append(" = ");
    error_msg.append(opt_values[0]);
    error_msg.append("\n\n");
    error_string.append(error_msg);
    return false;
  }
}

void Options::got_invalid_option(std::string &opt_name,
                                 std::string &error_string) {

  std::string error_msg;
  error_msg.append(opt_name);
  error_msg.append(": invalid option name");
  error_msg.append(". Check current options in config_template.cfg.");
  error_msg.append("\n\n");
  error_string.append(error_msg);
}

std::string Options::get_input_base_path(std::istream &is) {
  std::string line;
  std::string option;
  std::vector<std::string> values;
  std::string base_path;
  while (getline(is, line)) {
    if (line.find("INPUT_BASE_PATH") != std::string::npos) {
      if (TokenizeString(line, option, values)) // Not in a comment
      {
        base_path = values[0];
      }
    }
  }
  is.clear();
  is.seekg(0, std::ios_base::beg);
  return base_path;
}
