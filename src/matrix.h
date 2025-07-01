//
// Created by Даниил Маслов on 31.12.2022.
//

#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <iostream>
#include <array>
#include <iterator>
#include "biginteger.h"


template <size_t N, typename Field>
std::array<Field, N> &operator+=(std::array<Field, N> &that, const std::array<Field, N> &another) {
    for (size_t i = 0; i < N; ++i) {
        that[i] += another[i];
    }
    return that;
}

template <size_t N, typename Field>
std::array<Field, N> &operator-=(std::array<Field, N> &that, const std::array<Field, N> &another) {
    for (size_t i = 0; i < N; ++i) {
        that[i] -= another[i];
    }
    return that;
}

template <size_t N, typename Field>
std::array<Field, N> operator*=(std::array<Field, N> &that, const Field &lambda) {
    for (size_t i = 0; i < N; ++i) {
        that[i] *= lambda;
    }
    return that;
}

template <size_t N, typename Field>
std::array<Field, N> operator/=(std::array<Field, N> &that, const Field &lambda) {
    for (size_t i = 0; i < N; ++i) {
        that[i] /= lambda;
    }
    return that;
}

template <size_t N, typename Field>
std::array<Field, N> operator+(std::array<Field, N> that, const std::array<Field, N> &another) {
    that += another;
    return that;
}

template <size_t N, typename Field>
std::array<Field, N> operator-(std::array<Field, N> that, const std::array<Field, N> &another) {
    that -= another;
    return that;
}

template <size_t N, typename Field>
std::array<Field, N> operator*(std::array<Field, N> that, const Field &lambda) {
    that *= lambda;
    return that;
}

template <size_t N, typename Field>
std::array<Field, N> operator/(std::array<Field, N> that, const Field &lambda) {
    for (size_t i = 0; i < N; ++i) {
        that[i] /= lambda;
    }
    return that;
}

template <typename Field>
struct field_char {
    static const size_t value;
};

template <size_t N>
struct field_char<Residue<N>> {
    static const size_t value = N;
};

template <>
struct field_char<Rational> {
    static const size_t value = 0;
};

template <size_t N, size_t M, typename Field = Rational>
class Matrix {
  public:
    Matrix() {
        for (size_t i = 0; i < N; ++i) {
            matrix_[i].fill(Field());
        }
    }

    Matrix(const std::initializer_list<std::array<Field, M>> &another) {
        typename std::initializer_list<std::array<Field, M>>::iterator ptr = another.begin();
        for (size_t i = 0; ptr != another.end(); ++ptr, ++i) {
            matrix_[i] = *ptr;
        }
    }

    std::array<Field, M> &operator[](size_t index) {
        return matrix_[index];
    }

    const std::array<Field, M> &operator[](size_t index) const {
        return matrix_[index];
    }

    Matrix<N, M, Field> &operator+=(const Matrix<N, M, Field> &another) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                matrix_[i][j] += another.matrix_[i][j];
            }
        }
        return *this;
    }

    Matrix<N, M, Field> &operator-=(const Matrix<N, M, Field> &another) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                matrix_[i][j] -= another.matrix_[i][j];
            }
        }
        return *this;
    }

    Matrix<N, M, Field> &operator*=(const Field &lambda) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                matrix_[i][j] *= lambda;
            }
        };
        return *this;
    }

    Matrix<N, M, Field> &operator/=(const Field &lambda) {
        static_assert(is_field_v<Field>);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                matrix_[i][j] /= lambda;
            }
        }
        return *this;
    }

    Matrix<M, N, Field> transposed() const {
        Matrix<M, N, Field> ans;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                ans.matrix_[j][i] = matrix_[i][j];
            }
        }
        return ans;
    }

    Field det() const {
        static_assert(N == M);
        Matrix copy = *this;
        copy.triangulateMatrix();
        return copy.getDet();
    }

    size_t rank() const {
        Matrix copy = *this;
        return copy.triangulateMatrix();
    }

    Field trace() const {
        static_assert(N == M);
        return getTrace();
    }

    void invert() {
        static_assert(N == M);
        Matrix<N, 2 * N, Field> extendedMatrix;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                extendedMatrix.matrix_[i][j] = matrix_[i][j];
            }
        }
        for (size_t i = 0; i < N; ++i) {
            extendedMatrix.matrix_[i][N + i] = Field(1);
        }
        extendedMatrix.triangulateMatrix();
        for (size_t i = 1; i <= N; ++i) {
            const Field divider = extendedMatrix.matrix_[N - i][N - i];
            extendedMatrix.matrix_[N - i] /= divider;
            for (size_t j = 0; i + j < N; ++j) {
                extendedMatrix[j] -= extendedMatrix[N - i] * extendedMatrix[j][N - i];
            }
        }
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                matrix_[i][j] = extendedMatrix[i][N + j];
            }
        }
    }

    Matrix inverted() const {
        Matrix copy = *this;
        copy.invert();
        return copy;
    }

    std::array<Field, M> getRow(size_t index) const {
        return matrix_[index];
    }

    std::array<Field, N> getColumn(size_t index) const {
        std::array<Field, N> ans;
        for (size_t j = 0; j < N; ++j) {
            ans[j] = matrix_[j][index];
        }
        return ans;
    }

    bool operator==(const Matrix<N, M, Field> &another) const {
        bool ans = true;
        for (size_t i = 0; i < N && ans; ++i) {
            for (size_t j = 0; j < M && ans; ++j) {
                ans = (matrix_[i][j] == another.matrix_[i][j]);
            }
        }
        return ans;
    }

    void print() const {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                std::cerr << matrix_[i][j] << ' ';
            }
            std::cerr << '\n';
        }
    }

  protected:
    std::array<std::array<Field, M>, N> matrix_;

    template <size_t, size_t, typename> friend class Matrix;

    size_t triangulateMatrix() {
        static_assert(is_field_v<Field>);
        bool isChangedSign = false;
        size_t cnt = 0;
        for (size_t i = 0; i < M; ++i) {
            size_t j = cnt;
            while (j < N && matrix_[j][i] == Field()) {
                ++j;
            }
            if (j < N) {
                if (j != cnt) {
                    swap(matrix_[cnt], matrix_[j]);
                    isChangedSign = !isChangedSign;
                }
                j = 0;
                while (j < N) {
                    if (j != cnt) {
                        matrix_[j] -= matrix_[cnt] * (matrix_[j][i] / matrix_[cnt][i]);
                    }
                    ++j;
                }
                ++cnt;
            }
        }
        if (isChangedSign) {
            matrix_[N - 1] *= Field(-1);
        }
        return cnt;
    }

    Field getDet() const {
        static_assert(N == M);
        Field ans = 1;
        for (size_t i = 0; i < N; ++i) {
            ans *= matrix_[i][i];
        }
        return ans;
    }

    Field getTrace() const {
        Field ans = 0;
        for (size_t i = 0; i < N; ++i) {
            ans += matrix_[i][i];
        }
        return ans;
    }
};

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator+(Matrix<N, M, Field> a, const Matrix<N, M, Field> &b) {
    a += b;
    return a;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator-(Matrix<N, M, Field> a, const Matrix<N, M, Field> &b) {
    a -= b;
    return a;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(Matrix<N, M, Field> a, Field lambda) {
    a *= lambda;
    return a;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(Field lambda, Matrix<N, M, Field> a) {
    a *= lambda;
    return a;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator/(Matrix<N, M, Field> a, Field lambda) {
    a /= lambda;
    return a;
}

template <size_t A, size_t B, size_t C, size_t D, typename Field>
bool operator==(const Matrix<A, B, Field> &a, const Matrix<C, D, Field> &b) {
    if (A != C || B != D) {
        return false;
    }
    return a == b;
}

template <size_t N, typename Field>
Matrix<N, N, Field> &operator*=(Matrix<N, N, Field> &a, const Matrix<N, N, Field> &b) {
    Matrix<N, N, Field> ans;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            for (size_t k = 0; k < N; ++k) {
                ans[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    a = ans;
    return a;
}

template <size_t N, size_t M, size_t K, typename Field>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field> &a, const Matrix<M, K, Field> &b) {
    Matrix<N, K, Field> ans;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < K; ++j) {
            ans[i][j] = 0;
            for (size_t k = 0; k < M; ++k) {
                ans[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return ans;
}

template <size_t N, typename Field = Rational>
class SquareMatrix : public Matrix<N, N, Field> {
  public:
    SquareMatrix() : Matrix<N, N, Field>() {};

    SquareMatrix(const std::initializer_list<std::array<Field, N>> &another) : Matrix<N, N, Field>(another) {};
};

#endif //MATRIX_MATRIX_H