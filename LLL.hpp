#pragma once
#include "Basis.hpp"
#include "GramSchmidt.hpp"
#include <cmath>
#include <algorithm>

template <typename T, size_t dim, size_t emb>
bool is_size_reduced(const Basis<T, dim, emb>& b) {
    for (size_t i = 1; i < dim; ++i) {
        for (size_t j = 0; j < i; ++j) {
            T mu = b.GS_proj_coeffs[i][j];
            if (!(std::abs(mu) <= (0.5 + 1e-9))) {
                return false;
            }
        }
    }
    return true;
}

template <typename T, size_t dim, size_t emb>
bool lovasz_condition(const Basis<T, dim, emb>& b, T delta) {
    for (size_t k = 1; k < dim; ++k) {
        T norm_prev_sq = b.squared_norms[k - 1];
        T norm_curr_sq = b.squared_norms[k];
        T mu_k_prev = b.GS_proj_coeffs[k][k - 1];
        T lhs = delta * norm_prev_sq;
        T rhs = norm_curr_sq + (mu_k_prev * mu_k_prev * norm_prev_sq);
        if (!(lhs <= (rhs + 1e-9))) {
            return false;
        }
    }
    return true;
}

template <typename T, size_t dim, size_t emb>
bool is_lll_reduced(const Basis<T, dim, emb>& b, T delta) {
    return lovasz_condition(b, delta) && is_size_reduced(b);
}

template <typename T, size_t dim, size_t emb>
void update_gs_after_swap(Basis<T, dim, emb>& b, size_t k) {

    size_t km1 = k - 1;
 
    T mu = b.GS_proj_coeffs[k][km1];
    T B_k = b.squared_norms[k];
    T B_km1 = b.squared_norms[km1];

    T new_B_km1 = B_k + (mu * mu * B_km1);
    
    if (std::abs(new_B_km1) < 1e-9) return; 

    T mu_prime = mu * B_km1 / new_B_km1;
    
    b.squared_norms[km1] = new_B_km1;
    b.squared_norms[k] = (B_k * B_km1) / new_B_km1;
    b.GS_proj_coeffs[k][km1] = mu_prime;

    // Update the block of coefficients
    for (size_t j = 0; j < km1; ++j) {
        T temp = b.GS_proj_coeffs[km1][j];
        b.GS_proj_coeffs[km1][j] = b.GS_proj_coeffs[k][j];
        b.GS_proj_coeffs[k][j] = temp;
    }

    for (size_t i = k + 1; i < dim; ++i) {
        T t = b.GS_proj_coeffs[i][k];
        b.GS_proj_coeffs[i][k] = b.GS_proj_coeffs[i][km1] - mu * t;
        b.GS_proj_coeffs[i][km1] = t + mu_prime * b.GS_proj_coeffs[i][k];
    }
}

template <typename T, size_t dim, size_t emb>
Basis<T, dim, emb> lll(const Basis<T, dim, emb>& input_basis, T delta) {
   
    Basis<T, dim, emb> b = input_basis;

    gram_schmidt<false>(b); 

    size_t k = 1; 

    while (k < dim) {
        
        // --- Size Reduction ---
        for (size_t j = k - 1; (int) j >= 0; --j) {
            T mu_kj = b.GS_proj_coeffs[k][j];
            
            if (std::abs(mu_kj) > 0.5) {
                T qs = std::round(mu_kj);
                b[k] -= b[j] * qs;
                b.GS_proj_coeffs[k][j] -= qs;
                for (size_t l = 0; l < j; ++l) {
                    b.GS_proj_coeffs[k][l] -= qs * b.GS_proj_coeffs[j][l];
                }
            }
        }

        // --- Lovasz Condition Check ---
        T mu_k_km1 = b.GS_proj_coeffs[k][k - 1];
        T B_k = b.squared_norms[k];
        T B_km1 = b.squared_norms[k - 1];

        if (B_k >= (delta - mu_k_km1 * mu_k_km1) * B_km1) {
            k++;
        } else {
            b.swap(k, k - 1);
            update_gs_after_swap(b, k);
            k = std::max(k - 1, size_t(1));
        }
    }

    return b;
}

template <typename T, size_t dim, size_t emb>
T estimate_delta(const Basis<T, dim, emb>& b, T precision) {
    T best_delta = -1.0;
    
    for (T d = 0.25; d < 1.0; d += 0.01) {
        Basis<T, dim, emb> reduced = lll(b, d);
        if (is_lll_reduced(reduced, d)) {
            best_delta = d;
        } else {
            break; 
        }
    }
    
    if (best_delta > 0.0) {
        T start = best_delta;
        T end = std::min((T)1.0, start + 0.01);
        for(T d = start; d < end; d += precision) {
             Basis<T, dim, emb> reduced = lll(b, d);
             if (is_lll_reduced(reduced, d)) {
                 best_delta = d;
             }
        }
    }

    return best_delta;
}
