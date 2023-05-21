#include <iostream>
#include <vector>
#include <complex>
#include "matrix.hpp"

using namespace std;

int main()
{

	Matrix<3, 2> A = {
		73, 73, 
		69, 69,
		0, 1
	};

	Matrix<3, 1> B = {1, 1, 1};

	(A.kronecker_prod(B)).print();

	// system("pause");

	return 0;
}