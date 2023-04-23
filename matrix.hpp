#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>

using namespace std;

template<int m, int n, typename F = double> 
class Matrix {

private:
	F** data = new F*[m];

	template<class func> 
    void loop(func f); 

public:
	Matrix(); 
	Matrix(F** data);
	Matrix(vector<vector<F>> arr);

	inline F* operator [] (int i) const {return data[i];}
	Matrix operator + (Matrix const& other) const;
	Matrix operator - (Matrix const& other) const;
	Matrix operator * (F const& c) const;
	template<int l>
	Matrix<m,l,F> operator * (Matrix<n,l,F> const& other) const;

	void print();
};



//implementations

template<int m, int n, typename F> 
Matrix<m,n,F>::Matrix() {
		for (int i = 0; i < m; i++) {
			data[i] = new F[n];
		}
	}

template<int m, int n, typename F> 
Matrix<m,n,F>::Matrix(F** data) {
		*this = Matrix();

		auto set = [&] (int i, int j) {
			this->data[i][j] = data[i][j];
		};

		loop(set);
	}

template<int m, int n, typename F> 
Matrix<m,n,F>::Matrix(vector<vector<F>> arr) {
	if (arr.size() != m || arr[0].size() != n) {
		throw std::invalid_argument("Array dimensions must match matrix dimensions");
	}

	*this = Matrix();

	auto set = [&] (int i, int j) {
			this->data[i][j] = arr[i][j];
		};

	loop(set);

}

template<int m, int n, typename F>
template<int l>
Matrix<m,l,F> Matrix<m,n,F>::operator * (Matrix<n,l,F> const& other) const {
	F** new_data = new F*[m];
	for (int i = 0; i < m; i++) {
		new_data[i] = new F[l];
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < l; j++) {
			new_data[i][j] = 0;
			for (int k = 0; k < n; k++) {
				new_data[i][j] += data[i][k]*other[k][j];
			}
		}
	}

	return Matrix<m,l,F>(new_data);
}

template<int m, int n, typename F>
Matrix<m,n,F> Matrix<m,n,F>::operator + (Matrix<m,n,F> const& other) const {
	auto add = [&] (int i, int j) {
		this->data[i][j] += other[i][j];
	};

	loop(add);

	return *this;
}

template<int m, int n, typename F> 
template<class func>
void Matrix<m,n,F>::loop(func f) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				f(i,j);
			}
		}
	}

template<int m, int n, typename F> 
void Matrix<m,n,F>::print() {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				std::cout << data[i][j] << ' ';
			}
			std::cout << '\n';
		}
	}

#endif