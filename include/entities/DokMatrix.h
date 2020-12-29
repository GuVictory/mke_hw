//
// Created by GuVictory on 29.12.2020.
//

#ifndef DZ2_DOKMATRIX_H
#define DZ2_DOKMATRIX_H


#include "entities/AbstractMatrix.h"
#include "utils/Inaccuracy.h"

#include <map>

class DokMatrix : public AbstractMatrix {
public:
    using DataContainer = std::map<Index2d, Value>;

    explicit DokMatrix(Index side = 0);
    Index size() const override;
    Index2d shape() const override;
    Index non_zero() const noexcept;

    Value operator()(Index i, Index j) const;
    Value add(Index i, Index j, Value val);
    Value sub(Index i, Index j, Value val);
    Value set(Index i, Index j, Value val);
    const DataContainer& data() const noexcept;

private:
    template <class Callable>
    Value modify(Index2d pos, Callable&& f);

    DataContainer data_;
    Index m_;
};

template <class Callable>
Value DokMatrix::modify(Index2d ind, Callable&& f)
{
    Value result = 0;
    if (auto new_val = f(0); !iszero(new_val)) {
        auto [pos, inserted] = data_.insert({ind, new_val});
        if (!inserted) {
            result = pos->second;
            pos->second = f(pos->second);
        }
    }
    return result;
}


#endif //DZ2_DOKMATRIX_H
