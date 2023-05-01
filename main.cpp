#include <iostream>
#include <vector>
#include "matrix.hpp"

using namespace std;

int main() {
	
	vector<vector<double>> arr1 = {
		{1,0},
		{0,2},
		{0,0}
	};

	vector<vector<double>> arr2 = {
		{0,2,22},
		{0,0,1}
	};

	Matrix<3,2> A(arr1);
	Matrix<2,3> B(arr2);

	(A - A).print();

	return 0;
}