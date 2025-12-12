#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Basis.hpp"
#include "LLL.hpp"

template <typename T, size_t dim, size_t emb>
T calculate_volume(const Basis<T, dim, emb>& b) {
    T vol_sq = 1.0;
    for(size_t i = 0; i < dim; ++i) {
        vol_sq *= b.squared_norms[i];
    }
    return std::sqrt(vol_sq);
}

void print_status(const std::string& label, bool pass) {
    std::cout << label << ": " 
              << (pass ? "\033[32mPASS\033[0m" : "\033[31mFAIL\033[0m") 
              << "\n";
}

int main() {
    std::cout << "====================================\n";
    std::cout << "       LLL Implementation Test      \n";
    std::cout << "====================================\n\n";

    // --- Test Case Setup ---
    Basis<double, 3, 3> b;
    b[0] = {{1.0, 1.0, 1.0}};
    b[1] = {{-1.0, 0.0, 2.0}};
    b[2] = {{3.0, 5.0, 6.0}};
    double delta = 0.99;

    // --- Pre-Computation ---
    Basis<double, 3, 3> b_initial = b;
    gram_schmidt<false>(b_initial);
    double initial_vol = calculate_volume(b_initial);

    std::cout << "--- Initial Basis ---\n";
    b[0].print_vec("v0");
    b[1].print_vec("v1");
    b[2].print_vec("v2");
    std::cout << "Initial Volume: " << initial_vol << "\n\n";

    // --- Run LLL ---
    std::cout << "Running LLL(delta=" << delta << ")...\n";
    Basis<double, 3, 3> reduced = lll(b, delta);
    std::cout << "Done.\n\n";

    std::cout << "--- Reduced Basis ---\n";
    reduced[0].print_vec("v0");
    reduced[1].print_vec("v1");
    reduced[2].print_vec("v2");

    // --- Validation ---
    
    double final_vol = calculate_volume(reduced);
    double vol_diff = std::abs(initial_vol - final_vol);
    bool vol_pass = vol_diff < 1e-6;

    std::cout << "\n--- Verification Results ---\n";
    std::cout << "Final Volume:   " << final_vol << "\n";
    print_status("Volume Preserved", vol_pass);

    bool size_pass = is_size_reduced(reduced);
    print_status("Size Reduced    ", size_pass);

    bool lovasz_pass = lovasz_condition(reduced, delta);
    print_status("Lovasz Condition", lovasz_pass);
    
    double d = estimate_delta(b_initial, 0.00001);
    std::cout << "\n--- Estimated delta using heuristic ---\n";
    std::cout << "Delta:   " << d << "\n";

    return 0;
}
