#include <iostream>
#include <vector>
#include <complex>
#include "matrix.h"

using namespace std; 

int main()
{
	vector<vector<double>> arr = {
		{1, 0, 2},
		{2, 0, 1},
		{3, 1, 0}};

	const Matrix<3, 3> A(arr);

	Matrix<3, 3,std::complex<double>> B = {
		1.0 + 1i, 2         , 1,
		0       , 2.0 - 3.0i, 1i,
		1i      , 2.0       , 0};

	Matrix<3, 2> C = {
		1, 1,
		1, 1,
		1, 1};
	auto X = static_cast<Matrix<3,3,std::complex<double>>>(A);
	(B*Matrix<3,3,std::complex<double>>::id()).print();


	//system("pause");

	return 0;
}