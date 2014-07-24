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
void calc_source_vec(Matrix &nodes, IntMatrix &triangles, Matrix &source);
void SaveData(const char *filename, const int N, const double *X);
void SaveData(const char *filename, const int N, const int *X);
void get_mesh(int n, Matrix &nodes, IntMatrix &triangles);
void calc_stiff_matrix(Matrix &nodes, IntMatrix &triangles, Matrix &stiff_matrix);
void calc_bound_con(Matrix &nodes, IntMatrix &triangles, std::vector<int> &bnd);
}

using namespace easyfem;
int main()
{
	/*
	Matrix nodes(2, 3);
	IntMatrix T(10, 1);
	nodes(1, 1) = 0; nodes(2, 1) = 0;
	nodes(1, 2) = 2; nodes(2, 2) = 0;
	nodes(1, 3) = 2; nodes(2, 3) = 2;
	T(1) = 1; T(2) = 2; T(3) = 3;
	for(int i = 4; i <= 10; ++i)
	{
		T(i) = i - 3;
	}
	for(int i = 1; i <= 10; ++i)
	{
		printf("%d: %d\n", i, T(i, 1));
	}
	Matrix stiff_mat(0, 0);
	calc_stiff_matrix(nodes, T, stiff_mat);
	for(int i = 1; i <= stiff_mat.get_rows(); ++i)
	{
		for(int j = 1; j <= stiff_mat.get_cols(); ++j)
		{
			printf("%f,", stiff_mat(i, j));
		}
		printf("\n");
	}
	//*/
	/*
		Matrix nodes(2, 3);
		IntMatrix T(10,1);
		nodes(1, 1) = 0.; nodes(2, 1) = 0.;
		nodes(1, 2) = 1.; nodes(2, 2) = 0.;
		nodes(1, 3) = 0.; nodes(2, 3) = 1.;
		T(1) = 1; T(2) = 2; T(3) = 3;
		T(4) = 1; T(5) = 2; T(6) = 3; T(7) = 4;
		T(8) = 5; T(9) = 6; T(10) = 7;

		Matrix stiff_mat(0, 0);
		calc_stiff_matrix(nodes, T, stiff_mat);
		for(int i = 1; i <= stiff_mat.get_rows(); ++i)
		{
			for(int j = 1; j <= stiff_mat.get_cols(); ++j)
			{
				printf("%f,", stiff_mat(i, j));
			}
			printf("\n");
		}
	//*/
	/*
		Matrix nodes(0, 0);
		IntMatrix triangles(0, 0);
		get_mesh(2, nodes, triangles);


		for(int i = 1; i <= nodes.get_cols(); ++i)
		{
			printf("%d: (%f, %f)\n", i, nodes(1, i), nodes(2, i));
		}
		for(int i = 1; i <= triangles.get_cols(); ++i)
		{
			printf("%d: %d, %d, %d\n", i, triangles(1, i), triangles(2, i), triangles(3, i));
		}

		for(int i = 1; i <= triangles.get_cols(); ++i)
		{
			printf("%d:", i);
			for(int j = 4; j <= 10; ++j)
			{
				printf("%d,", triangles(j, i));
			}
			printf("\n");
		}
	//*/
	//*********************************************************
	//*
	Matrix nodes(0, 0);
	IntMatrix triangles(0, 0);
	Matrix stiff_matrix(0, 0);

	get_mesh(2, nodes, triangles);
	calc_stiff_matrix(nodes, triangles, stiff_matrix);
	int N = stiff_matrix.get_cols();
	SaveData("/home/lxy/workspace/A.txt", N * N, stiff_matrix.get_ptr(1));
	for(int i = 1; i <= 10; ++i)
	{
		printf("i: %d\n", triangles(i, 3));
	}
	/*
	Matrix source(0, 0);
	calc_source_vec(nodes, triangles, source);
	SaveData("/home/lbj-lxy/workspace/b.txt", source.size(), source.get_ptr(1));

	printf("here\n");
	//*/
	//*
	std::vector<int> bnd;
	calc_bound_con(nodes, triangles, bnd);
	printf("%ld\n", bnd.size());
	for(int i = 0; i < bnd.size(); ++i)
	{
		printf("%d: %d\n", i + 1, bnd[i]);
	}
	SaveData("/home/lxy/workspace/bnd.txt", bnd.size(), &bnd[0]);
	//*/
	printf("------------------Done!---------------------\n");
}
