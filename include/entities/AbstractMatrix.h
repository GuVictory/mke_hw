#ifndef DZ2_ABSTRACTMATRIX_H
#define DZ2_ABSTRACTMATRIX_H


#include "utils/Types.h"

#include <ostream>

class AbstractMatrix {
public:
    virtual ~AbstractMatrix();

    [[nodiscard]] virtual Index size() const;
    [[nodiscard]] virtual Index2d shape() const = 0;
};


#endif //DZ2_ABSTRACTMATRIX_H