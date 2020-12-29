#include "entities/CsrMatrix.h"
#include <utils/Inaccuracy.h>
#include <entities/Vec.h>

#include <algorithm>

using namespace std;

CsrMatrix::CsrMatrix()
        : data_()
        , indptr_(1)
        , indices_()
        , m_(0)
{
}

#define _(i, j) ((i) * (m_) + (j))

CsrMatrix::CsrMatrix(Index side, const DataContainer& vals)
        : diag_(side)
        , data_()
        , indptr_(side + 1)
        , indices_()
        , m_(side)
{
    if (!vals.empty()) {
        for (Index i = 0; i < m_; ++i) {
            for (Index j = 0; j < m_; ++j) {
                auto& v = vals[_(i, j)];
                if (!iszero(v)) {
                    if (i == j) {
                        diag_[i] = v;
                    } else {
                        data_.emplace_back(v);
                        indices_.emplace_back(j);
                        ++indptr_[i + 1];
                    }
                }
            }
        }

        data_.shrink_to_fit();
        indices_.shrink_to_fit();
        if (!data_.empty()) {
            for (Index i = 0; i < m_; ++i) {
                indptr_[i + 1] += indptr_[i];
            }
        }
    }
}

#undef _

CsrMatrix::CsrMatrix(const DokMatrix& dok)
        : CsrMatrix(dok.shape().first)
{
    data_.reserve(dok.non_zero());
    indices_.reserve(dok.non_zero());

    for (auto& p : dok.data()) {
        auto [ind, val] = p;
        if (!iszero(val)) {
            if (ind.first == ind.second) {
                diag_[ind.first] = val;
            } else {
                data_.emplace_back(val);
                indices_.emplace_back(ind.second);
                ++indptr_[ind.first + 1];
            }
        }
    }
    data_.shrink_to_fit();
    indices_.shrink_to_fit();

    if (!data_.empty()) {
        for (Index i = 0; i < m_; ++i) {
            indptr_[i + 1] += indptr_[i];
        }
    }
}

Index CsrMatrix::size() const
{
    return m_ * m_;
}

Index2d CsrMatrix::shape() const
{
    return {m_, m_};
}

Index CsrMatrix::non_zero() const
{
    return data_.size();
}


const CsrMatrix::DataContainer& CsrMatrix::diag() const noexcept {
    return diag_;
}

const CsrMatrix::DataContainer& CsrMatrix::data() const noexcept
{
    return data_;
}
const CsrMatrix::IndexContainer& CsrMatrix::indptr() const noexcept
{
    return indptr_;
}

const CsrMatrix::IndexContainer& CsrMatrix::indices() const noexcept
{
    return indices_;
}

Value CsrMatrix::operator()(Index i, Index j) const
{
    Value result = 0;

    if (i == j) {
        result = diag_[i];
    } else {
        auto b = begin(indices_);
        auto ifirst = b + indptr_[i];
        auto ilast = b + indptr_[i + 1];
        if (ifirst < ilast) {
            auto p = lower_bound(ifirst, ilast, j);
            if (p != ilast && *p == j) {
                result = data_[Index(p - b)];
            }
        }
    }

    return result;
}

Vec operator*(const CsrMatrix& lhs, const Vec& rhs)
{
    Vec result(lhs.m_, 0.);
    dot(result, lhs, rhs);
    return result;
}

Vec operator*(const Vec& lhs, const CsrMatrix& rhs)
{
    Vec result(rhs.m_, 0.);
    dot(result, lhs, rhs);
    return result;
}

void dot(Vec& result, const CsrMatrix& lhs, const Vec& rhs)
{
    for (Index i = 0; i < lhs.m_; ++i) {
        auto u = lhs.diag_[i] * rhs[i];
        auto k = lhs.indptr_[i];
        auto last = lhs.indptr_[i + 1];
        for (; k < last; ++k) {
            u += lhs.data_[k] * rhs[lhs.indices_[k]];
        }
        result[i] = u;
    }
}

void dot(Vec& result, const Vec& lhs, const CsrMatrix& rhs)
{
    for (Index i = 0; i < rhs.m_; ++i) {
        result[i] += lhs[i] * rhs.diag_[i];

        auto k = rhs.indptr_[i];
        auto last = rhs.indptr_[i + 1];
        for (; k < last; ++k) {
            auto j = rhs.indices_[k];
            result[j] += lhs[j] * rhs.data_[k];
        }
    }
}