#include <cmath>
#include <iostream>
#include "Basis.hpp"
#include "GramSchmidt.hpp"

int main() {
    Basis<double, 2, 3> b;

    b[0] = {{1.0, 1.0, 0.0}}; 
    b[1] = {{1.0, 0.0, 0.0}};

    std::cout << "--- Before Gram-Schmidt ---\n";
    b[0].print_vec();
    b[1].print_vec();
        
    gram_schmidt(b);

    std::cout << "\n--- After Gram-Schmidt ---\n";
    b[0].print_vec();
    b[1].print_vec();

    double dot_p = dot(b[0], b[1]);
    std::cout << "\nDot Product: " << dot_p << " (Should be 0.0)\n";
    std::cout << "Norm v0:     " << norm(b[0]) << " (Should be 1.0)\n";
    std::cout << "Norm v1:     " << norm(b[1]) << " (Should be 1.0)\n";

    return 0;
}
