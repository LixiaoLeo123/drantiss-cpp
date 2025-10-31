#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
namespace drantiss {
    class Matrix {
    private:
        static constexpr double SINGULAR_TOL = 1e-12;
        size_t m;
        size_t n;
        double* data;
        size_t index(size_t i, size_t j) const;
        Matrix inverse() const;
    public:
        Matrix();
        Matrix(size_t rows, size_t cols);
        Matrix(const Matrix& other);
        Matrix& operator=(const Matrix& other);
        ~Matrix() { delete[] data; };
        size_t rows() const { return m; };
        size_t cols() const { return n; };
        double& operator()(size_t i, size_t j);
        const double& operator()(size_t i, size_t j) const;
        Matrix operator+(const Matrix& other) const;
        Matrix operator-(const Matrix& other) const;
        Matrix operator*(const Matrix& other) const;
        Matrix operator*(double scalar) const;
        Matrix& operator+=(const Matrix& other);
        Matrix& operator-=(const Matrix& other);
        Matrix& operator*=(const Matrix& other);
        Matrix& operator*=(double scalar);
        Matrix operator~() const;
        double& operator[](size_t index);
        double operator[](size_t index) const;
        Matrix getMinor(size_t i, size_t j) const;
        double cofactor(size_t i, size_t j) const;
        double detByLaplace() const;
        double det() const;
        bool isSquare() const { return m == n; };
        bool isSingular() const;
        friend Matrix operator*(double scalar, const Matrix& mat);
        friend Matrix operator^(const Matrix& mat, int exponent);
        friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);
        friend std::istream& operator>>(std::istream& is, Matrix& mat);
        static Matrix I(size_t n);
        static Matrix diag(const double* values, size_t n);
        static Matrix colvec(const double* values, size_t n);
        void swapRows(size_t r1, size_t r2);
        void swapCols(size_t c1, size_t c2);
        void multiplyRow(size_t r, double k);
        void multiplyCol(size_t c, double k);
        void addRowMultiple(size_t r, size_t s, double k);
        void addColMultiple(size_t c, size_t d, double k);
    };
    inline size_t Matrix::index(size_t i, size_t j) const {
        return (i - 1) * n + j - 1;
    }
    inline Matrix::Matrix() : m(0), n(0), data(nullptr) {}
    inline Matrix::Matrix(size_t rows, size_t cols)
        : m(rows), n(cols), data(new double[rows * cols]()) {
    }
    inline bool Matrix::isSingular() const {
        return (abs(det()) < SINGULAR_TOL);
    }
    inline Matrix operator*(double scalar, const Matrix& mat) {
        return mat * scalar;
    }
}
#endif