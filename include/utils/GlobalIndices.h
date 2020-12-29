#ifndef DZ2_GLOBALINDICES_H
#define DZ2_GLOBALINDICES_H

#include "utils/Types.h"
#include "mke/Triangulation.h"

enum class Coord : Index { X, Y, Z };

template <class CoordIndex>
Index _g(const Index& node, const CoordIndex& coord, const Index& m) {
    return node * Triangulation::DIM + Index(coord);
}

template <class CoordIndex>
Index _v(const Index& node, const CoordIndex& coord) {
    return _g(node, coord, Triangulation::N);
}

template <class CoordIndex>
Index _s(const Index& node, const CoordIndex& coord) {
    return _g(node, coord, Triangulation::SN);
}

#endif //DZ2_GLOBALINDICES_H
