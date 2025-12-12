#pragma once
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>

template <typename T, size_t N> 
struct Vector {
    // T = type parameter, N = size parameter
    std::array<T, N> data;

    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }

    Vector<T, N> operator*(T scalar) const {
        Vector<T, N> result = *this;
        result *= scalar;
        return result;
    }

    Vector<T, N>& operator*=(T scalar) {
        for (size_t i = 0; i < N; ++i) {
            data[i] *= scalar;
        }
        return *this;
    }

    Vector<T, N> operator+(const Vector<T, N>& other) const {
        Vector<T, N> result = *this;
        result += other; 
        return result;
    }

    Vector<T, N>& operator+=(const Vector<T, N>& other) {
        for (size_t i = 0; i < N; ++i) {
            data[i] += other[i];
        }
        return *this;
    }

    Vector<T, N>& operator-=(const Vector<T, N>& other) {
        for (size_t i = 0; i < N; ++i) {
            data[i] -= other[i];
        }
        return *this;
    }

    Vector<T, N> operator-(const Vector<T, N>& other) const {
        Vector<T, N> result = *this;
        result -= other;
        return result;
    }

    void print_vec(const std::string& name = "", int precision = 4) {
        if (!name.empty()) std::cout << name << ": ";
        std::cout << "[ ";
        for (size_t i = 0; i < N; ++i) {
            std::cout 
                << std::fixed 
                << std::setprecision(precision) 
                << std::setw(8) 
                << data[i] << " ";
        }
        std::cout << "]\n";
    }

};

template <typename T, size_t N>
T dot(const Vector<T, N>& a, const Vector<T, N>& b) {
    T sum = 0;
    for (size_t i = 0; i < N; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

template <typename T, size_t N> 
Vector<T, N> operator*(T scalar, const Vector<T, N>& vec) {
    return vec * scalar;
}

template <typename T, size_t N>
long double norm(const Vector<T, N>& v) {
    return std::sqrt(dot(v,v));
}

template <typename T, size_t N>
long double norm_sq(const Vector<T, N>& v) {
    return dot(v,v);
}
