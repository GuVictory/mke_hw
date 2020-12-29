#ifndef DZ2_INACCURACY_H
#define DZ2_INACCURACY_H

#include "utils/Types.h"

enum class Tolerance { ZERO = 0, SINGLE = 1, DOUBLE = 2, TRIPLE = 3 };

bool isnear(Value lhs, Value rhs, Tolerance t = Tolerance::DOUBLE);

bool iszero(Value x, Tolerance t = Tolerance::DOUBLE);

Value sqr(Value x);

#endif //DZ2_INACCURACY_H
