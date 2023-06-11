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


int main()
{
	vec3 A = {0, 0, 1};

	vec3 B = {1,0,0};

	mat3 R = rotatexy<float>(PI/4);

	for (int i = 0; i < 9; i++) {
		std::cout << "B is: \n";
		mat3 X = mat3::id();
		(X).print();
		std::cout << "\n";
	}
	return 0;
}