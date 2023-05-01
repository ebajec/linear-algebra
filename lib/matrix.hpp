#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <initializer_list>
#include <vector>
#include "math.h"
#include <iostream>
#include <cassert>
#include "matrix.h"

using namespace std;
// constructors
template <int m, int n, typename F>
Matrix<m, n, F>::Matrix()
{
}

template <int m, int n, typename F>
Matrix<m, n, F>::Matrix(const initializer_list<F> &arr)
{
	int i = 0;
	for (F x : arr)
	{
		mem[i] = x;
		i++;
	}
}

template <int m, int n, typename F>
Matrix<m, n, F>::Matrix(F *new_data)
{
	*this = Matrix();

	for (int i = 0; i < m * n; i++)
	{
		this->mem[i] = new_data[i];
	}
}

template <int m, int n, typename F>
Matrix<m, n, F>::Matrix(const vector<vector<F>> &arr)
{
	if (arr.size() != m || arr[0].size() != n)
	{
		throw std::invalid_argument("Array dimensions must match matrix dimensions");
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			this->mem[i * n + j] = arr[i][j];
		}
	}
}

// operator overloads

template <int m, int n, typename F>
Matrix<m, n, F> Matrix<m, n, F>::operator+(const Matrix<m, n, F> &B) const
{
	Matrix<m, n, F> new_mat;
	for (int i = 0; i < m * n; i++)
	{
		new_mat.mem[i] += this->mem[i] + B.mem[i];
	}

	return new_mat;
}

template <int m, int n, typename F>
Matrix<m, n, F> Matrix<m, n, F>::operator-(const Matrix<m, n, F> &B) const
{
	Matrix<m, n, F> new_mat;
	for (int i = 0; i < m * n; i++)
	{
		new_mat.mem[i] += this->mem[i] - B.mem[i];
	}

	return new_mat;
}

template <int m, int n, typename F>
Matrix<m, n, F> Matrix<m, n, F>::operator*(F const &c)
{
	Matrix<m, n, F> new_mat;
	for (int i = 0; i < m * n; i++)
	{
		new_mat.mem[i] = mem[i] * c;
	}

	return new_mat;
}

template <int m, int n, typename F>
template <int l>
Matrix<m, l, F> Matrix<m, n, F>::operator*(const Matrix<n, l, F> &B) const
{
	F new_data[m * l] = {0};

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < l; j++)
		{
			for (int k = 0; k < n; k++)
			{
				new_data[i * l + j] += mem[i * n + k] * B[0][k * l + j];
			}
		}
	}
	return Matrix<m, l, F>(new_data);
}

template <int m, int n, typename F>
bool Matrix<m, n, F>::operator==(const Matrix<m, n, F> &other) const
{
	for (int i = 0; i < m * n; i++)
	{
		if (mem[i] != other[0][i])
		{
			return false;
		}
		return true;
	}
}

template <int m, int n, typename F>
bool Matrix<m, n, F>::operator!=(const Matrix<m, n, F> &other) const
{
	return !(*this == other);
}

// other
template <int m, int n, typename F>
void Matrix<m, n, F>::print()
{
	for (int i = 0; i < m * n; i++)
	{
		std::cout << mem[i] << ' ';

		if ((i + 1) % n == 0)
		{
			std::cout << '\n';
		}
	}
	std::cout << '\n';
}

template <int m, int n, typename F>
Matrix<1, n, F> Matrix<m, n, F>::row(int i) const
{
	Matrix<1, n, F> new_mat;

	for (int j = 0; j < n; j++)
	{
		new_mat[0][j] = mem[i * n + j];
	}

	return new_mat;
}

template <int m, int n, typename F>
Matrix<m, 1, F> Matrix<m, n, F>::col(int j) const
{
	Matrix<m, 1, F> new_mat;
	for (int i = 0; i < m; i++)
	{
		new_mat[0][i] = mem[i * n + j];
	}

	return new_mat;
}

template <int m, int n, typename F>
Matrix<m - 1, n - 1, F> Matrix<m, n, F>::submatrix(int row, int col) const
{
	assert((row < m) && (col < n));
	Matrix<m - 1, n - 1, F> new_mat;
	F *new_data = new_mat.data();

	int k = 0;

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (!(i == row || j == col))
			{
				new_data[k] = mem[i * n + j];
				k++;
			}
		}
	}

	return new_mat;
}

template <int m, int n, typename F>
Matrix<n, m, F> Matrix<m, n, F>::transpose() const
{
	Matrix<n, m, F> new_mat;

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			new_mat[0][j * m + i] = mem[i * n + j];
		}
	}

	return new_mat;
}

/******NON MEMBER FUNCTIONS*****/

/*Compute determinant via laplacian expansion*/
template <int n, typename F>
F det_laplace(Matrix<n, n, F> A)
{
	if (n == 1)
	{
		return A[0][0];
	}

	F value = 0;
	for (int i = 0; i < n; i++)
	{
		value += pow(-1, i) * (A[0][i] * det_laplace(A.submatrix(0, i)));
	}

	return value;
}

#endif