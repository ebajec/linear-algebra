#pragma once
#ifndef MATRIX_HPP
#define MATRIX_HPP
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

    inline F *operator[](int i) { return mem + i * n; }
    inline const F *operator[](int i) const { return mem + i * n; }
    Matrix<m, n, F> operator+(const Matrix<m, n, F> &B) const;
    Matrix<m, n, F> operator-(const Matrix<m, n, F> &B) const;
    Matrix<m, n, F> operator*(const F &c) const;
    template <int l>
    Matrix<m, l, F> operator*(const Matrix<n, l, F> &B) const;
    Matrix<m, n, F> operator^(const int &pow) const;
    /*Creates a new matrix by placing the two inputs side by side.*/
    template <int l>
    Matrix<m, n + l, F> operator|(const Matrix<m, l, F> &other) const;
    bool operator==(const Matrix<m, n, F> &other) const;
    bool operator!=(const Matrix<m, n, F> &other) const;

    template <typename T>
    operator Matrix<m, n, T>() const
    {
        T new_data[m * n];
        for (int i = 0; i < m * n; i++)
        {
            new_data[i] = static_cast<T>(mem[i]);
        }
        return Matrix<m, n, T>(new_data);
    }

    F *data() { return &(this->mem[0]); }
    void print() const;
    int rows() { return m; }
    int cols() { return n; }

    Matrix<1, n, F> row(int i) const;
    Matrix<m, 1, F> col(int j) const;
    Matrix<m - 1, n - 1, F> submatrix(int i, int j) const;
    Matrix<n, m, F> transpose() const;
    template <int k, int l>
    Matrix<m + k, n + l, F> direct_sum(const Matrix<k, l, F> &other) const;
    template <int k, int l>
    Matrix<m*k,n*l,F> kronecker_prod(const Matrix<k,l,F> &other) const;

    static Matrix<m, n, F> id();

    // Don't really know how to avoid implementing this outside of class.
    // This is so you can multiply from the right as well.
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

// needed to do this for base cases if matrix size gets reduced in recursive
// function call
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

/*****NON MEMBER STUFF*****/

template <int n, typename F>
F det_laplace(Matrix<n, n, F> A);

#include "matrix.tpp"

#endif