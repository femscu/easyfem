/*
 * main.cpp
 *
 *  Created on: Mar 18, 2014
 *      Author: lbj-lxy
 */

#include <vector>
#include <cstdio>
#include <cmath>

#include "../INCLUDE/matrix.h"

using namespace easyfem;
namespace easyfem
{
void nablaw(Matrix &xx);
void SaveData(const char *filename, const int N, const double *X);
void SaveData(const char *filename, const int N, const int *X);
void get_mesh(int n, Matrix &P, IntMatrix &T, IntMatrix &E);
void calc_stiff_matrix(Matrix &nodes, IntMatrix &triangles, Matrix &stiff_matrix);
void read_mesh(Matrix &P, IntMatrix &T);
}

using namespace easyfem;
void test4()
{
	Matrix P(2, 3);
	P(1, 1) =  0.;  P(2, 1) = 0.;
	P(1, 2) =  1.;  P(2, 2) = 0.2;
	P(1, 3) = .5;    P(2, 3) = 0.5;
	IntMatrix T(10, 1);
	T(1) = 1; T(2) = 2; T(3) = 3;
	for(int i = 4; i <= 10; ++i)
	{
		T(i) = i - 3;
	}

	Matrix A(0, 0);
	calc_stiff_matrix(P, T, A);
	for(int i = 1; i <= A.get_num_rows(); ++i)
	{
		for(int j = 1; j <= A.get_num_cols(); ++j)
		{
			printf("%10f", A(i,j));
		}
		printf("\n");
	}
}
void test3()
{
	Matrix P(0, 0);
	IntMatrix T(0, 0);

	read_mesh(P, T);
	for(int j = 1; j <= P.get_num_cols(); ++j)
	{
		printf("%d: %f, %f\n", j, P(1, j), P(2, j));
	}

	puts("------------------------------------");
	for(int j = 1; j <= T.get_num_cols(); ++j)
	{
		printf("%d: ", j);
		for(int i = 1; i <= 10; ++i)
		{
			printf("%5d", T(i, j));
		}
		printf("\n");
	}

	return;
}
void test2()
{
	Matrix P(0, 0);
	IntMatrix T(0, 0);
	Matrix stiff_matrix(0, 0);

	read_mesh(P, T);
	puts("-----------------P------------");
	for(int j = 1; j <= P.get_num_cols(); ++j)
	{
		printf("%d: %f, %f\n", j, P(1, j), P(2, j));
	}
	puts("-----------------T-----------------");
	for(int j = 1; j <= T.get_num_cols(); ++j)
	{
		printf("%d: ", j);
		for(int i = 1; i <= 10; ++i)
		{
			printf("%5d", T(i, j));
		}
		printf("\n");
	}
	//*/
	SaveData("/home/ljm/workspace/lbj.txt", 10 * T.get_num_cols(), T.get_ptr(1));
	calc_stiff_matrix(P, T, stiff_matrix);
	int dof = stiff_matrix.get_num_cols();
	SaveData("/home/ljm/workspace/stiff.txt", dof * dof, stiff_matrix.get_ptr(1));
}

int main()
{

	test4();
	printf("------------------Done!---------------------\n");
}
