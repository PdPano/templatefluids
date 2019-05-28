#ifndef TIME_INTEGRATOR_TOOL_FACTORY_HPP
#define TIME_INTEGRATOR_TOOL_FACTORY_HPP

#include <memory>

class TimeIntegratorTool;
class Options;
struct PointFunctions;
class CartesianGrid;

std::shared_ptr<TimeIntegratorTool> create_time_integrator_tool(
    Options& opt, PointFunctions& pf, CartesianGrid& grid);

#endif /* TIME_INTEGRATOR_TOOL_FACTORY_HPP */
