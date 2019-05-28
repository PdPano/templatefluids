#include "nan_checker.hpp"

double default_if_nan(double v, double def)
{
    if (v != v) {
        return def;
    }
    return v;
}
