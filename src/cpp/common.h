#pragma once
#ifndef CHALC_COMMON_H
#define CHALC_COMMON_H

#include <vector>
#include <memory>
#include <tuple>
#include <map>
#include <algorithm>
#include <numeric>

namespace chalc {
    typedef double value_t;
    typedef long long int index_t;
    namespace common {
        using std::vector, std::map, std::tuple, std::iota, std::tie, std::shared_ptr, std::make_shared;
    }
}

#endif