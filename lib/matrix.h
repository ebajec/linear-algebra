#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <initializer_list>
#include <vector>
#include <type_traits>

using namespace std;
template <int m, int n, typename F = double>
class Matrix
{
protected:
    F mem[m * n] = {0};

public:
    Matrix();
    Matrix(const Matrix<m, n, F> &);
    Matrix(F *new_data);
    Matrix(const vector<vector<F>> &arr);
    Matrix(const initializer_list<F> &arr);

    inline F *operator[](int i) { return mem + i * m; }
    inline const F *operator[](int i) const { return mem + i * m; }
    Matrix<m, n, F> operator+(const Matrix<m, n, F> &B) const;
    Matrix<m, n, F> operator-(const Matrix<m, n, F> &B) const;
    Matrix<m, n, F> operator*(const F &c);
    template <int l>
    Matrix<m, l, F> operator*(const Matrix<n, l, F> &B) const;
    Matrix<m, n, F> operator^(const int &pow) const;
    bool operator==(const Matrix<m, n, F> &other) const;
    bool operator!=(const Matrix<m, n, F> &other) const;

    F *data() { return &(this->mem[0]); }
    void print() const;
    int rows() { return m; }
    int cols() { return n; }

    Matrix<1, n, F> row(int i) const;
    Matrix<m, 1, F> col(int j) const;
    Matrix<m - 1, n - 1, F> submatrix(int i, int j) const;

    Matrix<n, m, F> transpose() const;

    static Matrix<m,n,F> id()
    {
        static_assert(m == n, "Matrix must be square");
        Matrix<m, n, F> identity;
        for (int i = 0; i < n; i++)
        {
            identity[0][i * n + i] = 1;
        }
        return identity;
    }

    // Don't really know how to avoid implementing this outside of class
    friend Matrix<m, n, F> operator*(const F &c, const Matrix<m, n, F> &B)
    {
        F new_data[m * n];

        for (int i = 0; i < m * n; i++)
        {
            new_data[i] = B.mem[i] * c;
        }

        return Matrix<m, n>(new_data);
    }
};

template <int n, typename F>
class Matrix<0, n, F>
{
public:
    F *operator[](int i) const
    {
        return nullptr;
    }
    F *data() { return nullptr; }
    Matrix<0, n, F> submatrix(int i, int j) const { return Matrix<0, n>(); }
};

template <int n, typename F>
class Matrix<n, 0, F>
{
public:
    F *operator[](int i) const
    {
        return nullptr;
    }
    F *data() { return nullptr; }
    Matrix<n, 0, F> submatrix(int i, int j) const { return Matrix<n, 0>(); }
};

template <typename F>
class Matrix<0, 0, F>
{
public:
    F *operator[](int i) const
    {
        return nullptr;
    }
    F *data() { return nullptr; }
    Matrix<0, 0, F> submatrix(int i, int j) const { return Matrix<0, 0>(); }
};

// template <int m, int n, typename F = double, int >
// class MatrixGenerator

/*****NON MEMBER STUFF*****/

template <int n, typename F>
F det_laplace(Matrix<n, n, F> A);

#include "matrix.hpp"

#endif