#pragma once
#include "Basis.hpp"
#include <cmath>

template <bool Normalize = true, typename T, size_t dim, size_t emb>
void gram_schmidt(Basis<T, dim, emb>& b) {
    
    for (size_t k = 0; k < dim; ++k) b.GS_proj_coeffs[k][k] = 1.0;

    for (size_t i = 0; i < dim; ++i) {
        
        T norm_sq = dot(b[i], b[i]);

        if constexpr (Normalize) {
            T len = std::sqrt(norm_sq);
            b[i] *= (T(1.0) / len); 
            b.squared_norms[i] = 1.0; 
        } else {
            b.squared_norms[i] = norm_sq;
        }

        for (size_t j = i + 1; j < dim; ++j) {
            T proj_coeff = dot(b[j], b[i]);
            
            if constexpr (!Normalize) { proj_coeff /= b.squared_norms[i]; }
            
            b.GS_proj_coeffs[j][i] = proj_coeff;
            b[j] -= b[i] * proj_coeff;
        }
    }
}
