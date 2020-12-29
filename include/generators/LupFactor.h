#ifndef DZ2_LUPFACTOR_H
#define DZ2_LUPFACTOR_H


#include "entities/DenseMatrix.h"
#include "entities/Vec.h"
#include <cmath>

class LupFactor {
public:
    using PivotType = std::vector<Index>;

    explicit LupFactor(const DenseMatrix& in);
    explicit LupFactor(DenseMatrix&& in);

    const DenseMatrix& mat() const;
    const PivotType& pivot() const;

    LupFactor& factor();
    Vec solve(const Vec& v) const;
    DenseMatrix inverse() const;
    Value det() const;

private:
    Index m_;
    DenseMatrix mat_;
    PivotType pivot_;
};


#endif //DZ2_LUPFACTOR_H
