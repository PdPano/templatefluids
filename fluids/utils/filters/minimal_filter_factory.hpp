#ifndef MINIMAL_FILTER_FACTORY_HPP
#define MINIMAL_FILTER_FACTORY_HPP

#include <memory>
class MinimalFilter;
std::shared_ptr<MinimalFilter> create_minimal_filter(int order);

#endif /* MINIMAL_FILTER_FACTORY_HPP */
