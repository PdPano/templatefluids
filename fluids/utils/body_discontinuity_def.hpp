#ifndef BODY_DISCONTINUITY_DEF_HPP
#define BODY_DISCONTINUITY_DEF_HPP
#include "generic_discontinuity_def.hpp"

class BodyDiscontinuity : public GenericDiscontinuity {
public:
    BodyDiscontinuity()
        : GenericDiscontinuity({}, {})
        , T(1.0)
    {
    }
    double T;
};

#endif /* BODY_DISCONTINUITY_DEF_HPP */
