#include "entities/AbstractMatrix.h"

using namespace std;

AbstractMatrix::~AbstractMatrix() = default;

Index AbstractMatrix::size() const
{
    auto s = shape();
    return s.first * s.second;
}
