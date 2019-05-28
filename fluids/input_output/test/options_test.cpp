#include "../options.hpp"
#include "gtest/gtest.h"

std::string sample_file = ("INPUT_MESHFILE_NAME = input/gridInfo.dat\n"
                           "INPUT_IMMERSED_INTERFACE = immersedInterface2.dat\n"
                           "INPUT_FLOW_CONFIGURATION = initialCondition.dat\n"
                           "OUTPUT_FILE_NAME = result\n"
                           "OUTPUT_BASE_PATH = output/\n"
                           "INPUT_TYPE = DEFAULT\n"
                           "OUTPUT_TYPE = VTK_OUT\n"
                           " \n"
                           "GAMMA = 1.4\n"
                           "REYNOLDS = 1\n"
                           "MACH = 0.4\n"
                           "PRANDTL = 1\n"
                           "\n"
                           "SOLVER_TYPE = KARAGIOZIS\n"
                           "INTEGRATOR_TYPE = RUNGE_KUTTA_SEMI_IMPLICIT\n"
                           "MIX_PARAM = 0.1\n"
                           "CFL = 0.2\n"
                           "T_MAX=1.0\n"
                           "PRINT_INTERVAL=0.05\n");

class OptionsTest : public testing::Test {
protected:
    virtual void SetUp()
    {
        std::istringstream stream(sample_file);
        opts = std::make_unique<Options>();
        opts_parse = std::make_unique<Options>(stream);
    }
    std::unique_ptr<Options> opts;
    std::unique_ptr<Options> opts_parse;
};

TEST_F(OptionsTest, testDefaultConstructor) { ASSERT_EQ(opts->gam(), 1.4); }

TEST_F(OptionsTest, testFileParsing)
{
    ASSERT_EQ(opts_parse->reynolds(), 1.0);
    ASSERT_EQ(opts_parse->input_flow_configuration_file_name(),
        "./input/initialCondition.dat");
}
