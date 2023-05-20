#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <initializer_list>
#include <vector>
#include "math.h"
#include <iostream>
#include <cassert>
#include "matrix.h"

#define DBL_EPSILON 2.2204460492503131e-016

using namespace std;

// constructors
template <int m, int n, typename F>
Matrix<m, n, F>::Matrix()
{
}

template <int m, int n, typename F>
Matrix<m, n, F>::Matrix(const Matrix<m, n, F> &other)
{
	for (int i = 0; i < m * n; i++)
	{
		mem[i] = other[0][i];
	}
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
Matrix<m, n, F> Matrix<m, n, F>::operator*(F const &c) const
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
	// need to do this because Matrix<n, l, F> is a different type and access
	// to private fields is not allowed.
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
template <int l>
inline Matrix<m, n + l, F> Matrix<m, n, F>::operator|(const Matrix<m, l, F> &other) const
{
	F new_data[m * (l + n)] = {0};

	const int cols = n + l;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			new_data[i * cols + j] = mem[i * n + j];
		}
		for (int j = 0; j < l; j++)
		{
			new_data[i * cols + j + n] = other[0][i * l + j];
		}
	}

	return Matrix<m, n + l, F>(new_data);
}

template <int m, int n, typename F>
Matrix<m, n, F> Matrix<m, n, F>::operator^(const int &pow) const
{
	static_assert(m == n, "Matrix must be square");
	if (pow < 0)
	{
		throw std::invalid_argument("Negative matrix powers are not always defined");
	}

	Matrix<m, n, F> new_mat = Matrix<m, n, F>::id();

	for (int i = 0; i < pow; i++)
	{
		new_mat = new_mat * (*this);
	}

	return new_mat;
}

template <int m, int n, typename F>
bool Matrix<m, n, F>::operator==(const Matrix<m, n, F> &other) const
{
	for (int i = 0; i < m * n; i++)
	{
		if (mem[i] != other[0][i])

			return false;
	}
	return true;
}

template <int m, int n, typename F>
bool Matrix<m, n, F>::operator!=(const Matrix<m, n, F> &other) const
{
	return !(*this == other);
}

// other
template <int m, int n, typename F>
void Matrix<m, n, F>::print() const
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

	if (m == 1 || n == 1)
	{
		return new_mat;
	}

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

template <int m, int n, typename F>
template <int k, int l>
Matrix<m + k, n + l, F> Matrix<m, n, F>::direct_sum(const Matrix<k, l, F> &other) const
{
	F new_data[(m + k) * (n + l)] = {0};

	// This looks gross, but I would have to use the modulo operator a lot
	// to do it differently. I feel like this would be faster.
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			new_data[i * (n + l) + j] = mem[i * n + j];
		}
	}
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < l; j++)
		{
			new_data[(i + m) * (n + l) + j + n] = other[0][i * l + j];
		}
	}

	return Matrix<m + k, n + l, F>(new_data);
}

template <int m, int n, typename F>
Matrix<m, n, F> Matrix<m, n, F>::id()
{
	static_assert(m == n, "Matrix must be square");

	Matrix<m, n, F> identity;
	for (int i = 0; i < n; i++)
	{
		identity[0][i * n + i] = 1;
	}
	return identity;
}

/******NON MEMBER FUNCTIONS*****/

template <int m, int n, typename F>
Matrix<m, n, F> gauss_elim(Matrix<m, n, F> A)
{
	F data[m * n];
	for (int i = 0; i < m * n; i++)
	{
		data[i] = A[0][i];
	}

	// add c*r1 to r2
	auto row_add = [&data](int r1, int r2, F c)
	{
		for (int j = 0; j < n; j++)
		{
			data[r2 * n + j] += c * data[r1 * n + j];
		}
		return;
	};
	// swap r1 and r2
	auto row_swap = [&data](int r1, int r2)
	{
		for (int i = 0; i < n; i++)
		{
			F temp = data[r1 * n + i];
			data[r1 * n + i] = data[r2 * n + i];
			data[r2 * n + i] = temp;
		}
		return;
	};

	int row_nonzero = 0;
	for (int col = 0; col < n; col++)
	{
		F *top_entry = data + row_nonzero * n + col;

		// this attempts to swap rows and make top_entry nonzero
		if (abs(*top_entry) <= DBL_EPSILON)
		{
			for (int row = row_nonzero; row < m; row++)
			{
				if (abs(data[row * n + col]) > DBL_EPSILON)
				{
					row_swap(row_nonzero, row);
					break;
				}
			}
		}

		// we only continue if last part was successful
		F top_entry_val = *top_entry;
		if (abs(top_entry_val) > DBL_EPSILON)
		{
			// now we add row_nonzero to the ones below to get zeroes
			for (int row = row_nonzero + 1; row < m; row++)
			{
				F val = data[row * n + col];
				if (abs(val) > DBL_EPSILON)
				{
					row_add(row_nonzero, row, -val * pow(top_entry_val, -1));
				}
			}
			// go to next row once we're done
			row_nonzero++;
		}
	}

	return Matrix<m, n, F>(&data[0]);
}

/*Compute determinant via laplacian expansion*/
template <int n, typename F>
F det_laplace(Matrix<n, n, F> A)
{
	if (n == 0)
	{
		return 1;
	}
	if (n == 1)
	{
		return A[0][0];
	}
	else
	{
		F value = 0;

		for (int i = 0; i < n; i++)
		{
			value += pow(-1, i) * (A[0][i] * det_laplace(A.submatrix(0, i)));
		}
		return value;
	}
}

template <int n, typename F>
Matrix<n, n, F> inv(Matrix<n, n, F> A)
{
	F det = det_laplace(A);

	if (det > DBL_EPSILON)
	{
	}

	return pow(det, -1) * adj(A);
}

/*
 * Computes adjucate matrix.
 */
template <int n, typename F>
Matrix<n, n, F> adj(Matrix<n, n, F> A)
{

	if (n == 1)
	{
		return Matrix<n, n, F>::id();
	}

	Matrix<n, n, F> adjucate;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			adjucate[j][i] = pow(-1, i + j) * det_laplace(A.submatrix(i, j));
		}
	}
	return adjucate;
}

#endif