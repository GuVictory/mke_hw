#ifndef DZ2_LOCALEQGEN_H
#define DZ2_LOCALEQGEN_H

#include <entities/Vec.h>
#include <entities/DenseMatrix.h>
#include <mke/Triangulation.h>

#include <functional>

using LocalEqGen = std::function<std::pair<DenseMatrix, Vec>(
        const Triangulation&, const Triangulation::FiniteElement&)>;

LocalEqGen::result_type gen_local(const Triangulation& t,
                                  const Triangulation::FiniteElement& elem);


#endif //DZ2_LOCALEQGEN_H
