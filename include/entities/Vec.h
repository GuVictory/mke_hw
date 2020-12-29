#ifndef DZ2_VEC_H
#define DZ2_VEC_H


#include <ostream>
#include <vector>
#include "utils/Types.h"

class Vec : public std::vector<Value> {
    using Super = std::vector<Value>;

public:
    using Super::Super;

    // Перегрузим основные операции
    Vec& operator+=(const Vec& rhs);
    Vec& operator-=(const Vec& rhs);
    Vec& operator*=(Value rhs);
    Vec& operator/=(Value rhs);

    friend std::ostream& operator<<(std::ostream& os, const Vec& obj);

    friend bool operator==(const Vec& lhs, const Vec& rhs);
    friend bool operator!=(const Vec& lhs, const Vec& rhs);

    friend Vec operator+(const Vec& lhs, const Vec& rhs);
    friend Vec operator-(const Vec& lhs, const Vec& rhs);
    friend Vec operator*(const Vec& lhs, Value rhs);
    friend Vec operator*(Value lhs, const Vec& rhs);
    friend Vec operator/(const Vec& lhs, Value rhs);

    friend Value dot(const Vec& lhs, const Vec& rhs);
    friend Value cnorm(const Vec& lhs);
    friend Value cdist(const Vec& lhs, const Vec& rhs);
    friend Value sqr(const Vec& obj);

};

#endif //DZ2_VEC_H
