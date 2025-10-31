#include "Matrix.h"
#include "Essentials.h"
namespace drantiss {
    Matrix::Matrix(const Matrix& other)
        : m(other.m), n(other.n), data(new double[other.m * other.n]) {
        size_t total = m * n;
        for (size_t i = 0; i < total; ++i) {
            data[i] = other.data[i];
        }
    }
    Matrix& Matrix::operator=(const Matrix& other) {
        if (this != &other) {
            delete[] data;
            m = other.m;
            n = other.n;
            data = new double[m * n];
            size_t total = m * n;
            for (size_t i = 0; i < total; ++i) {
                data[i] = other.data[i];
            }
        }
        return *this;
    }
    double& Matrix::operator()(size_t i, size_t j) {
        if (i > m || j > n || i == 0 || j == 0) {
            std::cerr << "Error: Matrix index (" << i << ", " << j << ") out of range [1, "
                << m << "] x [1, " << n << "]\n";
            std::exit(EXIT_FAILURE);
        }
        return data[index(i, j)];
    }
    const double& Matrix::operator()(size_t i, size_t j) const {
        if (i > m || j > n || i == 0 || j == 0) {
            std::cerr << "Error: Matrix index (" << i << ", " << j << ") out of range [1, "
                << m << "] x [1, " << n << "]\n";
            std::exit(EXIT_FAILURE);
        }
        return data[index(i, j)];
    }

    Matrix Matrix::operator+(const Matrix& other) const {
        if (m != other.m || n != other.n) {
            std::cerr << "Error: Matrix dimensions mismatch for addition: ("
                << m << "x" << n << ") vs (" << other.m << "x" << other.n << ")\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix result(m, n);
        size_t total = m * n;
        for (size_t i = 0; i < total; ++i) {
            result.data[i] = data[i] + other.data[i];
        }
        return result;
    }

    Matrix Matrix::operator-(const Matrix& other) const {
        if (m != other.m || n != other.n) {
            std::cerr << "Error: Matrix dimensions mismatch for subtraction: ("
                << m << "x" << n << ") vs (" << other.m << "x" << other.n << ")\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix result(m, n);
        size_t total = m * n;
        for (size_t i = 0; i < total; ++i) {
            result.data[i] = data[i] - other.data[i];
        }
        return result;
    }

    Matrix Matrix::operator*(const Matrix& other) const {
        if (n != other.m) {
            std::cerr << "Error: Matrix multiplication dimension mismatch: ("
                << m << "x" << n << ") * (" << other.m << "x" << other.n << ")\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix result(m, other.n);
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < other.n; ++j) {
                double sum = 0.0;
                for (size_t r = 0; r < n; ++r) {
                    sum += (*this)(i + 1, r + 1) * other(r + 1, j + 1);
                }
                result(i + 1, j + 1) = sum;
            }
        }
        return result;
    }
    Matrix Matrix::operator*(double scalar) const {
        Matrix result(m, n);
        size_t total = m * n;
        for (size_t i = 0; i < total; ++i) {
            result.data[i] = data[i] * scalar;
        }
        return result;
    }
    Matrix& Matrix::operator+=(const Matrix& other) {
        if (m != other.m || n != other.n) {
            std::cerr << "Error: Matrix dimensions mismatch for +=\n";
            std::exit(EXIT_FAILURE);
        }
        size_t total = m * n;
        for (size_t i = 0; i < total; ++i) {
            data[i] += other.data[i];
        }
        return *this;
    }

    Matrix& Matrix::operator-=(const Matrix& other) {
        if (m != other.m || n != other.n) {
            std::cerr << "Error: Matrix dimensions mismatch for -=\n";
            std::exit(EXIT_FAILURE);
        }
        size_t total = m * n;
        for (size_t i = 0; i < total; ++i) {
            data[i] -= other.data[i];
        }
        return *this;
    }
    Matrix& Matrix::operator*=(const Matrix& other) {
        if (n != other.m) {
            std::cerr << "Error: Matrix multiplication dimension mismatch: ("
                << m << "x" << n << ") * (" << other.m << "x" << other.n << ")\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix result = (*this) * other;
        delete[] data;
        m = result.m;
        n = result.n;
        data = result.data;
        result.data = nullptr;
        return *this;
    }
    Matrix& Matrix::operator*=(double scalar) {
        size_t total = m * n;
        for (size_t i = 0; i < total; ++i) {
            data[i] *= scalar;
        }
        return *this;
    }
    Matrix Matrix::operator~() const {
        Matrix result(n, m);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                result.data[j * m + i] = data[i * n + j];
        return result;
    }
    double& Matrix::operator[](size_t index) {
        if (index == 0 || index > m * n) {
            std::cerr << "Error: Linear index " << index
                << " out of range [1, " << m * n << "].\n";
            std::exit(EXIT_FAILURE);
        }
        return data[index - 1];
    }

    double Matrix::operator[](size_t index) const {
        if (index < 1 || index > m * n) {
            std::cerr << "Error: Linear index " << index
                << " out of range [1, " << m * n << "].\n";
            std::exit(EXIT_FAILURE);
        }
        return data[index - 1];
    }
    Matrix Matrix::getMinor(size_t i, size_t j) const {
        if (i > m || j > n || i == 0 || j == 0) {
            std::cerr << "Error: Minor index (" << i << ", " << j << ") out of range [1, "
                << m << "] x [1, " << n << "].\n";
            std::exit(EXIT_FAILURE);
        }
        size_t new_m = m - 1;
        size_t new_n = n - 1;
        Matrix result(new_m, new_n);
        if (new_m == 0 || new_n == 0) {
            return result;
        }
        size_t r = 0;
        for (size_t ii = 0; ii < m; ++ii) {
            if (ii == i - 1) continue;
            size_t c = 0;
            for (size_t jj = 0; jj < n; ++jj) {
                if (jj == j - 1) continue;
                result.data[r * new_n + c] = data[ii * n + jj];
                ++c;
            }
            ++r;
        }
        return result;
    }
    double Matrix::cofactor(size_t i, size_t j) const {
        if (m != n) {
            std::cerr << "Error: Cofactor only defined for square matrices.\n";
            std::exit(EXIT_FAILURE);
        }
        if (i > m || j > n || i == 0 || j == 0) {
            std::cerr << "Error: Cofactor index (" << i << ", " << j << ") out of range [1, "
                << m << "] x [1, " << n << "].\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix minor = getMinor(i, j);
        double det = minor.det();
        return ((i + j) % 2 == 0 ? 1 : -1) * det;
    }
    double Matrix::detByLaplace() const {
        if (m != n) {
            std::cerr << "Error: Determinant only defined for square matrices.\n";
            std::exit(EXIT_FAILURE);
        }
        if (m == 0) return 1.0;
        if (m == 1) return data[0];
        if (m == 2) return data[0] * data[3] - data[1] * data[2];
        double result = 0.0;
        for (size_t j = 1; j <= n; ++j) {
            Matrix minor = getMinor(1, j);
            double cofactor = ((j % 2 == 1) ? 1.0 : -1.0) * minor.det();
            result += data[j - 1] * cofactor;
        }
        return result;
    }
    double Matrix::det() const {
        if (m != n) {
            std::cerr << "Error: Determinant only defined for square matrices.\n";
            std::exit(EXIT_FAILURE);
        }
        if (m == 0) return 1.0;
        if (m == 1) return data[0];
        double* A = new double[n * n];
        size_t total = n * n;
        for (size_t i = 0; i < total; ++i) {
            A[i] = data[i];
        }
        int swaps = 0;
        for (size_t k = 0; k < n; ++k) {
            size_t pivot_row = k;
            double max_val = abs(A[k * n + k]);
            for (size_t i = k + 1; i < n; ++i) {
                double abs_val = abs(A[i * n + k]);
                if (abs_val > max_val) {
                    max_val = abs_val;
                    pivot_row = i;
                }
            }
            if (max_val < SINGULAR_TOL) {
                delete[] A;
                return 0.0;
            }
            if (pivot_row != k) {
                for (size_t j = 0; j < n; ++j) {
                    swap(A[k * n + j], A[pivot_row * n + j]);
                }
                swaps++;
            }
            for (size_t i = k + 1; i < n; ++i) {
                double factor = A[i * n + k] / A[k * n + k];
                for (size_t j = k + 1; j < n; ++j) {
                    A[i * n + j] -= factor * A[k * n + j];
                }
            }
        }
        double detU = 1.0;
        for (size_t i = 0; i < n; ++i) {
            detU *= A[i * n + i];
        }
        double sign = (swaps % 2 == 0) ? 1.0 : -1.0;
        delete[] A;
        return sign * detU;
    }
    Matrix operator^(const Matrix& mat, int exponent) {
        if (mat.m != mat.n) {
            std::cerr << "Error: Matrix power only defined for square matrices.\n";
            std::exit(EXIT_FAILURE);
        }
        size_t n = mat.m;
        if (exponent == 0) {
            Matrix I(n, n);
            for (size_t i = 0; i < n; ++i)
                I.data[i * n + i] = 1.0;
            return I;
        }
        if (exponent < 0) {
            Matrix inv = mat.inverse();
            return inv ^ (-exponent);
        }
        Matrix result = Matrix(n, n);
        for (size_t i = 0; i < n; ++i)
            result.data[i * n + i] = 1.0;
        Matrix base = mat;
        int exp = exponent;
        while (exp > 0) {
            if (exp & 1)
                result *= base;
            base = base * base;
            exp >>= 1;
        }
        return result;
    }
    Matrix Matrix::inverse() const {
        if (m != n) {
            std::cerr << "Error: Inverse only for square matrices.\n";
            std::exit(EXIT_FAILURE);
        }
        double d = det();
        if (abs(d) < SINGULAR_TOL) {
            std::cerr << "Error: Matrix is singular (det = 0), no inverse.\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix adj(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                adj.data[j * n + i] = cofactor(i + 1, j + 1);
            }
        }
        return adj * (1.0 / d);
    }
    
    std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
        for (size_t i = 0; i < mat.m; ++i) {
            for (size_t j = 0; j < mat.n; ++j) {
                os << mat.data[mat.index(i + 1, j + 1)];
                if (j < mat.n - 1) {
                    os << " ";
                }
            }
            os << "\n";
        }
        return os;
    }
    std::istream& operator>>(std::istream& is, Matrix& mat) {
        size_t total = mat.m * mat.n;
        for (size_t i = 0; i < total; ++i) {
            is >> mat.data[i];
        }
        return is;
    }
    Matrix Matrix::I(size_t n) {
        Matrix mat(n, n);
        for (size_t i = 0; i < n; ++i) {
            mat.data[i * n + i] = 1;
        }
        return mat;
    }
    Matrix Matrix::diag(const double* values, size_t n) {
        if (values == nullptr && n > 0) {
            std::cerr << "Error: diag() received null pointer with n = " << n << "\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix mat(n, n);
        for (size_t i = 0; i < n; ++i) {
            mat.data[i * n + i] = values[i];
        }
        return mat;
    }
    Matrix Matrix::colvec(const double* values, size_t n) {
        if (values == nullptr && n > 0) {
            std::cerr << "Error: colvec() received null pointer with n = " << n << "\n";
            std::exit(EXIT_FAILURE);
        }
        Matrix mat(n, 1);
        for (size_t i = 0; i < n; ++i) {
            mat.data[i] = values[i];
        }
        return mat;
    }
    void Matrix::swapRows(size_t r1, size_t r2) {
        if (r1 > m || r2 > m || r1 == 0 || r2 == 0) {
            std::cerr << "Error: swapRows index " << r1 << " or " << r2 << " out of range [1, "
                << m << "].\n";
            return;
        }
        for (size_t j = 1; j <= n; ++j)
            swap(data[index(r1, j)], data[index(r2, j)]);
    }
    void Matrix::swapCols(size_t c1, size_t c2) {
        if (c1 > n || c2 > n || c1 == 0 || c2 == 0) {
            std::cerr << "Error: swapCols index " << c1 << " or " << c2 << " out of range [1, "
                << n << "].\n";
            return;
        }
        for (size_t j = 1; j <= m; ++j)
            swap(data[index(j, c1)], data[index(j, c2)]);
    }
    void Matrix::multiplyRow(size_t r, double k) {
        if (r > m || r == 0) {
            std::cerr << "Error: multiplyRow index " << r << " out of range [1, "
                << m << "].\n";
            return;
        }
        for (size_t j = 1; j <= n; ++j)
            data[index(r, j)] *= k;
    }
    void Matrix::multiplyCol(size_t c, double k) {
        if (c > m || c == 0) {
            std::cerr << "Error: multiplyCow index " << c << " out of range [1, "
                << m << "].\n";
            return;
        }
        for (size_t j = 1; j <= m; ++j)
            data[index(j, c)] *= k;
    }
    void Matrix::addRowMultiple(size_t r, size_t s, double k) {
        if (r > m || s > m || r == 0 || s == 0) {
            std::cerr << "Error: addRowMultiple index " << r << " or " << s << " out of range [1, "
                << m << "].\n";
            return;
        }
        for (size_t j = 1; j <= n; ++j)
            data[index(r, j)] += k * data[index(s, j)];
    }
    void Matrix::addColMultiple(size_t c, size_t d, double k) {
        if (c > m || d > m || c == 0 || d == 0) {
            std::cerr << "Error: addColMultiple index " << c << " or " << d << " out of range [1, "
                << n << "].\n";
            return;
        }
        for (size_t j = 1; j <= n; ++j)
            data[index(j, c)] += k * data[index(j, d)];
    }

}
