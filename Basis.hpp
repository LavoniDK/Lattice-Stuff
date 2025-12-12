#pragma once
#include <cmath>
#include "Vector.hpp"

template <typename T, size_t dim, size_t emb>
struct Basis {
    std::array<Vector<T, emb>, dim> vectors;
    T GS_proj_coeffs[dim][dim]; // Gram-Schmidt projection coefficients
    T squared_norms[dim]; // Squared norms

    Vector<T, emb>& operator[](size_t i) { return vectors[i]; }
    const Vector<T, emb>& operator[](size_t i) const { return vectors[i]; }

    void swap(size_t i, size_t j) {
        std::swap(vectors[i], vectors[j]);
    }
};
