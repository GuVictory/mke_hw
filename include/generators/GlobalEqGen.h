#ifndef DZ2_GLOBALEQGEN_H
#define DZ2_GLOBALEQGEN_H

#include <entities/CsrMatrix.h>
#include "generators/LocalEqGen.h"
#include "mke/Triangulation.h"

std::pair<DenseMatrix, Vec> GlobalEqGenDense(const Triangulation& t,
                                                      LocalEqGen gen);


std::pair<CsrMatrix, Vec> GlobalEqGenSystem(const Triangulation& t,
                                              LocalEqGen gen);

#endif //DZ2_GLOBALEQGEN_H
