#pragma once

#include <algorithm>
#include <vector>

template <size_t N, size_t M, typename T = int64_t>
class Matrix {
 public:
  Matrix(std::vector<std::vector<T>> vec) {
    std::vector v_temp(M, T());
    elem_.resize(N);
    for (size_t i = 0; i < N; ++i) {
      elem_[i] = v_temp;
    }
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        (*this)(i, j) = vec[i][j];
      }
    }
  }
  Matrix(T elem) {
    std::vector v_temp(M, T());
    elem_.resize(N);
    for (size_t i = 0; i < N; ++i) {
      elem_[i] = v_temp;
    }
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        (*this)(i, j) = elem;
      }
    }
  }
  Matrix() {
    std::vector v_temp(M, T());
    elem_.resize(N);
    for (size_t i = 0; i < N; ++i) {
      elem_[i] = v_temp;
    }
  }
  bool operator==(const Matrix<N, M, T>& right) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        if (elem_[i][j] != right(i, j)) {
          return false;
        }
      }
    }
    return true;
  }

  Matrix<M, N, T> Transposed() {
    Matrix<M, N, T> res;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        res(i, j) = elem_[j][i];
      }
    }
    return res;
  }
  const T& operator()(size_t row, size_t col) const { return elem_[row][col]; }
  T& operator()(size_t row, size_t col) { return elem_[row][col]; }
  Matrix<N, M, T>& operator+=(const Matrix<N, M, T>& right) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        (*this)(i, j) += right(i, j);
      }
    }
    return *this;
  }
  Matrix<N, M, T> operator+(const Matrix<N, M, T>& right) {
    Matrix res(*this);
    res += right;
    return res;
  }
  Matrix<N, M, T>& operator-=(const Matrix<N, M, T>& right) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        (*this)(i, j) -= right(i, j);
      }
    }
    return *this;
  }
  Matrix<N, M, T> operator-(const Matrix<N, M, T>& right) {
    Matrix res(*this);
    res -= right;
    return res;
  }
  Matrix<N, M, T> operator*(T num) {
    Matrix<N, M, T> res;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        res(i, j) += elem_[i][j] * num;
      }
    }
    return res;
  }

 private:
  std::vector<std::vector<T>> elem_;
};
template <size_t N, size_t M, size_t K, typename T = int64_t>
Matrix<N, K, T> operator*(const Matrix<N, M, T> kLeft,
                          const Matrix<M, K, T> kRight) {
  Matrix<N, K, T> res;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t step = 0; step < M; ++step) {
        res(i, j) += kLeft(i, step) * kRight(step, j);
      }
    }
  }
  return res;
}
template <size_t N, typename T>
class Matrix<N, N, T> {
 public:
  Matrix() {
    std::vector<T> v_temp(N, T());
    elem_.resize(N);
    for (size_t i = 0; i < N; ++i) {
      elem_[i] = v_temp;
    }
  }
  Matrix(const std::vector<std::vector<T>>& vec) {
    std::vector<T> v_temp(N, T());
    elem_.resize(N);
    for (size_t i = 0; i < N; ++i) {
      elem_[i] = v_temp;
    }
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        elem_[i][j] = vec[i][j];
      }
    }
  }
  Matrix(T elem) {
    std::vector<T> v_temp(N, T());
    elem_.resize(N);
    for (size_t i = 0; i < N; ++i) {
      elem_[i] = v_temp;
    }
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        elem_[i][j] = elem;
      }
    }
  }
  Matrix<N, N, T> operator*(const T kLeft) {
    Matrix<N, N, T> res;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        res(i, j) = elem_[i][j] * kLeft;
      }
    }
    return res;
  }
  Matrix<N, N, T>& operator+=(const Matrix<N, N, T>& right) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        elem_[i][j] += right(i, j);
      }
    }
    return *this;
  }
  Matrix<N, N, T> operator+(const Matrix<N, N, T>& right) {
    Matrix res(*this);
    res += right;
    return res;
  }
  Matrix<N, N, T>& operator-=(const Matrix<N, N, T>& right) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        elem_[i][j] -= right(i, j);
      }
    }
    return *this;
  }
  Matrix<N, N, T> operator-(const Matrix<N, N, T>& right) {
    Matrix res(*this);
    res -= right;
    return res;
  }
  bool operator==(const Matrix<N, N, T>& right) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        if (elem_[i][j] != right(i, j)) {
          return false;
        }
      }
    }
    return true;
  }
  T& operator()(size_t row, size_t col) { return elem_[row][col]; }
  const T& operator()(size_t row, size_t col) const { return elem_[row][col]; }
  Matrix<N, N, T> Transposed() {
    Matrix res;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        res(i, j) = elem_[j][i];
      }
    }
    return res;
  }
  T Trace() {
    T res = 0;
    for (size_t i = 0; i < N; ++i) {
      res += elem_[i][i];
    }
    return res;
  }

 private:
  std::vector<std::vector<T>> elem_;
};
