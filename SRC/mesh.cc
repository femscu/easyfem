/*
 * mesh.cc
 *
 *  Created on: Mar 18, 2014
 *      Author: lbj-lxy
 */
#include <cstdio>
#include "../INCLUDE/matrix.h"

namespace easyfem
{
void get_mesh(int n, Matrix &nodes, IntMatrix &triangles)
{
	Matrix x(n + 1, 1);
	Matrix y(n + 1, 1);

	for(int i = 1; i <= n + 1; ++i)
	{
		x(i) = -1 + (i - 1) * 2. / n;
		y(i) = x(i);
	}

	nodes.resize(2, (n + 1) * (n + 1));
	for(int j = 1, count = 1; j <= n + 1; ++j)
	{
		for(int i = 1; i <= n + 1; ++i)
		{
			nodes(2, count) = y(j);
			nodes(1, count) = x(i);
			++count;
		}
	}

	triangles.resize(10, 2 * n * n);
	for(int j = 1, count = 1; j <= n; ++j)
	{
		for(int i = 1; i <= n; ++i)
		{
			int index1 = (j - 1) * (n + 1) + i;
			int index2 = index1 + 1;
			int index3 = index2 + n + 1;
			int index4 = index3 - 1;
			triangles(1, count) = index1; triangles(2, count) = index2; triangles(3, count) = index3;
			++count;
			triangles(1, count) = index1; triangles(2, count) = index3; triangles(3, count) = index4;
			++count;
		}
	}

	int dof = 1;
	for(int i = 1; i <= triangles.get_cols(); ++i)
	{
		triangles(4, i) = dof;
		++dof;
	}

	for(int i = 1; i <= triangles.get_cols(); ++i)
	{
		triangles(5, i) = dof;
		triangles(6, i) = dof + 1;
		triangles(7, i) = dof + 2;
		triangles(8, i) = dof + 3;
		triangles(9, i) = dof + 4;
		triangles(10, i) = dof + 5;
		dof += 6;
	}
}
}//namespace easyfem


