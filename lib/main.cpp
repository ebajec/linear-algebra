#include <iostream>
#include <vector>
#include "matrix.h"

using namespace std;

int main()
{
	vector<vector<double>> arr = {
		{1, 0, 2},
		{2, 0, 1},
		{3, 1, 0}};

	const Matrix<3, 3> A(arr);

	Matrix<2, 3> B = {
		0, 2, 1,
		0, 1, 0};

	Matrix<3, 2> C = {
		1, 1,
		1, 1,
		1, 1};

	C.print();
	B.print();

	auto X = C * B;

	X.print();
	B.transpose().print();
	std::cout << det_laplace(X);

	//system("pause");

	return 0;
}