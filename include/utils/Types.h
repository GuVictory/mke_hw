#ifndef DZ2_TYPES_H
#define DZ2_TYPES_H

#include <cstdlib>
#include <ostream>
#include <utility>

using Value = double;
using Index = size_t;
using Index2d = std::pair<Index, Index>;

std::ostream& operator<<(std::ostream& os, const Index2d& obj);

#endif //DZ2_TYPES_H
