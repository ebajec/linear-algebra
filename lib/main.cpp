/*This file is just for running tests

I might add some better tests
*/

#include <iostream>
#include <vector>
#include <complex>
#include "matrix.hpp"

#define PI 3.14159265359

using namespace std;

using vec3 = matrix<3, 1>;
using mat3 = matrix<3, 3>;

vec3 n = vec3::random(10);
vec3 l = vec3::random(10);
vec3 l_o = vec3::random(10);

template <int m, int n, typename F>
bool test_rref(matrix<m,n,F> A) {
	std::cout << "\nRREF test running:\n";
	std::cout << "original is:\n";
	A.print();
	std::cout << "RREF is:\n";
	auto RREF = rref(A);
	RREF.print();
	std::cout << "solution is:\n";
	auto sol = RREF.col(n-1);
	sol.print();

	std::cout << "difference is\n";
	auto diff = A.col(n-1) - matrix<m,m,F>(A)*sol;
	diff.print();

	bool success = dot(diff,diff) < 1.192092896e-07F;
	return success;
}

int main()
{
	vec3 A = {0, 0, 1};

	matrix<3,4> B = {
		1,0,1,1,
		3,5,0,0,
		4,0,2,3
		};

	test_rref(B);



	return 0;
}