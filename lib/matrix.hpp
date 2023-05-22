#pragma once
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <initializer_list>
#include <vector>
#include <type_traits>

using namespace std;

template <int m, int n, typename F = double>
class matrix
{
protected:
    F mem[m * n] = {0};

public:
    matrix();
    matrix(const matrix<m, n, F> &);
    matrix(F *new_data);
    matrix(const vector<vector<F>> &arr);
    matrix(const initializer_list<F> &arr);

    inline F *operator[](int i) { return mem + i * n; }
    inline const F *operator[](int i) const { return mem + i * n; }
    matrix<m, n, F> operator+(const matrix<m, n, F> &B) const;
    matrix<m, n, F> operator-(const matrix<m, n, F> &B) const;
    matrix<m, n, F> operator*(const F &c) const;
    template <int l>
    matrix<m, l, F> operator*(const matrix<n, l, F> &B) const;
    matrix<m, n, F> operator^(const int &pow) const;
    /*Creates a new matrix by placing the two inputs side by side.*/
    template <int l>
    matrix<m, n + l, F> operator|(const matrix<m, l, F> &other) const;
    bool operator==(const matrix<m, n, F> &other) const;
    bool operator!=(const matrix<m, n, F> &other) const;

    template <typename T>
    operator matrix<m, n, T>() const
    {
        T new_data[m * n];
        for (int i = 0; i < m * n; i++)
        {
            new_data[i] = static_cast<T>(mem[i]);
        }
        return matrix<m, n, T>(new_data);
    }

    F *data() { return &(this->mem[0]); }
    void print() const;
    int rows() { return m; }
    int cols() { return n; }

    matrix<1, n, F> row(int i) const;
    matrix<m, 1, F> col(int j) const;
    matrix<m - 1, n - 1, F> submatrix(int i, int j) const;
    matrix<n, m, F> transpose() const;
    template <int k, int l>
    matrix<m + k, n + l, F> direct_sum(const matrix<k, l, F> &other) const;
    template <int k, int l>
    matrix<m*k,n*l,F> kronecker_prod(const matrix<k,l,F> &other) const;

    static matrix<m, n, F> id();

    // Don't really know how to avoid implementing this outside of class.
    // This is so you can multiply from the right as well.
    friend matrix<m, n, F> operator*(const F &c, const matrix<m, n, F> &B)
    {
        F new_data[m * n];

        for (int i = 0; i < m * n; i++)
        {
            new_data[i] = B.mem[i] * c;
        }

        return matrix<m, n>(new_data);
    }
};

// needed to do this for base cases if matrix size gets reduced in recursive
// function call
template <int n, typename F>
class matrix<0, n, F>
{
public:
    F *operator[](int i) const
    {
        return nullptr;
    }
    F *data() { return nullptr; }
    matrix<0, n, F> submatrix(int i, int j) const { return matrix<0, n>(); }
};

template <int n, typename F>
class matrix<n, 0, F>
{
public:
    F *operator[](int i) const
    {
        return nullptr;
    }
    F *data() { return nullptr; }
    matrix<n, 0, F> submatrix(int i, int j) const { return matrix<n, 0>(); }
};

template <typename F>
class matrix<0, 0, F>
{
public:
    F *operator[](int i) const
    {
        return nullptr;
    }
    F *data() { return nullptr; }
    matrix<0, 0, F> submatrix(int i, int j) const { return matrix<0, 0>(); }
};

/*****NON MEMBER STUFF*****/

template <int n, typename F>
F det_laplace(matrix<n, n, F> A);

#include "matrix.tpp"

#endif