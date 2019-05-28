#include "writer.hpp"

Writer::Writer(Options& opt)
    : base_name(opt.output_file_name())
    , output_type(opt.output_type())
    , counter(opt.output_counter())
    , pf(PointFunctions(opt.mach(), opt.gam()))
{
}
