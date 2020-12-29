#ifndef DZ2_CSRMATRIX_H
#define DZ2_CSRMATRIX_H


#include <array>
#include <iostream>
#include <iterator>
#include <set>
#include <type_traits>
#include <vector>

#include "entities/AbstractMatrix.h"
#include "entities/DokMatrix.h"
#include "entities/Vec.h"
#include "utils/Inaccuracy.h"

class CsrMatrix : public AbstractMatrix {
public:
    using DataContainer = std::vector<Value>;
    using IndexContainer = std::vector<Index>;

    CsrMatrix();
    explicit CsrMatrix(Index side, const DataContainer& vals = DataContainer());
    explicit CsrMatrix(const DokMatrix& dok);

    [[nodiscard]] Index size() const override;
    [[nodiscard]] Index2d shape() const override;
    [[nodiscard]] Index non_zero() const;

    [[nodiscard]] const DataContainer& diag() const noexcept;
    [[nodiscard]] const DataContainer& data() const noexcept;
    [[nodiscard]] const IndexContainer& indptr() const noexcept;
    [[nodiscard]] const IndexContainer& indices() const noexcept;

    Value operator()(Index i, Index j) const;
    friend Vec operator*(const CsrMatrix& lhs, const Vec& rhs);
    friend Vec operator*(const Vec& lhs, const CsrMatrix& rhs);
    friend void dot(Vec& result, const CsrMatrix& lhs, const Vec& rhs);
    friend void dot(Vec& result, const Vec& lhs, const CsrMatrix& rhs);

    friend Vec solve(const CsrMatrix& lhs, const Vec& rhs, Vec x0);
    friend Vec psolve(const CsrMatrix& lhs, const Vec& rhs, Vec x0);

private:
    DataContainer diag_;
    DataContainer data_;
    IndexContainer indptr_;
    IndexContainer indices_;
    Index m_;
};


#endif //DZ2_CSRMATRIX_H
