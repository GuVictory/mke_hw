#include "generators/LocalEqGen.h"

#include "utils/Constants.h"
#include <utils/Inaccuracy.h>
#include "utils/GlobalIndices.h"

#include <generators/LupFactor.h>

#include <iterator>

using namespace std;

class V17InternalGen {
    const Triangulation::FiniteElement::Data data;

public:
    explicit V17InternalGen(const Triangulation::FiniteElement::Data& data);
    [[nodiscard]] DenseMatrix gk() const;
    [[nodiscard]] static DenseMatrix internal(Value val) ;

private:
    [[nodiscard]] DenseMatrix dq() const;
};

class V17SurfaceGen {
    const Triangulation::SurfaceElement::Data data;

public:
    explicit V17SurfaceGen(const Triangulation::SurfaceElement::Data& data);

    [[nodiscard]] DenseMatrix boundary(Value val) const;
};

class V17Gen {
    const Triangulation& t;
    const Triangulation::FiniteElement elem;
    const V17InternalGen igen;

public:
    V17Gen(const Triangulation& t, const Triangulation::FiniteElement& elem);

    [[nodiscard]] DenseMatrix gk() const;
    [[nodiscard]] DenseMatrix internal(Value val) const;
    [[nodiscard]] DenseMatrix boundary(Value val) const;
    [[nodiscard]] static DenseMatrix boundary(const Triangulation::SurfaceElement& face, Index i,
                         Value val) ;
};

//! Генератор, который позволяет совместить матрицу по границу и объему
LocalEqGen::result_type gen_local(const Triangulation& t,
                                  const Triangulation::FiniteElement& elem)
{
    V17Gen gen(t, elem);

    auto vk = elem.volume();
    auto e3 = DenseMatrix::eye(Triangulation::DIM);
    auto s0 = kroneker_product(gen.internal(vk * consts::rho / 20), e3);
    auto sb = kroneker_product(gen.boundary(1. / 24), e3);

    Vec tk(Triangulation::DIM * Triangulation::N, 0.);
    for (Index i = 0; i < Triangulation::N; ++i) {
        tk[_v(i, Coord::Z)] = -consts::p;
    }

    auto mat = vk * gen.gk() - sqr(consts::omega) * s0 + consts::mu * sb;
    auto vec = sb * tk;

    return {mat, vec};
}


V17InternalGen::V17InternalGen(const Triangulation::FiniteElement::Data& data_)
        : data(data_)
{
}

DenseMatrix V17InternalGen::dq() const
{
    DenseMatrix mat({Triangulation::N, Triangulation::N}, 1);
    for (Index i = 0; i < Triangulation::DIM; ++i) {
        for (Index j = 0; j < Triangulation::N; ++j) {
            mat(i, j) = data[j][i];
        }
    }
    auto r = LupFactor(mat).factor();

    array<Vec, Triangulation::DIM> f = {
            r.solve({1, 0, 0, 0}), r.solve({0, 1, 0, 0}), r.solve({0, 0, 1, 0})};

    DenseMatrix result(
            {2 * Triangulation::DIM, Triangulation::N * Triangulation::DIM});
    for (Index i = 0; i < Triangulation::N; ++i) {
        auto j = Triangulation::DIM * i;
        result(0, j + 0) = f[0][i];
        result(1, j + 1) = f[1][i];
        result(2, j + 2) = f[2][i];
        result(3, j + 0) = f[1][i];
        result(3, j + 1) = f[0][i];
        result(4, j + 0) = f[2][i];
        result(4, j + 2) = f[0][i];
        result(5, j + 1) = f[2][i];
        result(5, j + 2) = f[1][i];
    }

    return result;
}

DenseMatrix V17InternalGen::gk() const
{
    static DenseMatrix cm({2 * Triangulation::DIM, 2 * Triangulation::DIM},
                          begin(consts::C), end(consts::C));

    auto m = dq();
    auto mt = m.transpose();
    return mt * cm * m;
}

DenseMatrix V17InternalGen::internal(Value val)
{
    DenseMatrix result({Triangulation::N, Triangulation::N}, val);
    for (Index i = 0; i < Triangulation::N; ++i) {
        result(i, i) *= 2;
    }
    return result;
}

V17SurfaceGen::V17SurfaceGen(const Triangulation::SurfaceElement::Data& data_)
        : data(data_)
{
}

DenseMatrix V17SurfaceGen::boundary(Value val) const
{
    auto s = norm(cross(data[0] - data[2], data[1] - data[2]));

    DenseMatrix result({Triangulation::SN, Triangulation::SN}, val * s);
    for (Index i = 0; i < Triangulation::SN; ++i) {
        result(i, i) *= 2;
    }
    return result;
}

V17Gen::V17Gen(const Triangulation& t_,
               const Triangulation::FiniteElement& elem_)
        : t(t_)
        , elem(elem_)
        , igen(elem_.data())
{
}
DenseMatrix V17Gen::internal(Value val) const
{
    return igen.internal(val);
}

DenseMatrix V17Gen::boundary(Value val) const
{
    DenseMatrix result({4, 4});
    for (Index i = 0; i < Triangulation::N; ++i) {
        auto f = elem.face(i);
        if (t.on_third(f)) {
            result += boundary(f, i, val);
        }
    }
    return result;
}

DenseMatrix V17Gen::boundary(const Triangulation::SurfaceElement& face,
                             Index i, Value val)
{
    DenseMatrix result({4, 4}, 0.);
    auto m = (i + 1) % Triangulation::N;
    auto n = (i + 2) % Triangulation::N;
    auto p = (i + 3) % Triangulation::N;
    auto c = 2 * face.area();

    result(m, m) = result(n, n) = result(p, p) = 2 * c * val;
    result(m, n) = result(m, p) = result(n, p) = result(n, m) = result(p, m)
            = result(p, n) = 1 * c * val;

    return result;
}

DenseMatrix V17Gen::gk() const
{
    return igen.gk();
}
