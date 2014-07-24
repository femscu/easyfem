/*
 * mesh.cc
 *
 *  Created on: Mar 18, 2014
 *      Author: lbj-lxy
 */
#include <cstdio>
#include <vector>
#include "../INCLUDE/matrix.h"

namespace easyfem
{

int find(int N, IntMatrix &E, int node1, int node2)
{
	for(int i = 1; i <= N; ++i)
	{
		if((E(1, i) == node1) && (E(2, i) == node2))
		{
			return i;
		}

		if((E(1, i) == node2) && (E(2, i) == node1))
		{
			return i;
		}
	}

	return 0;
}


// P contains coordinates of nodals of the mesh
// T contains the information of all triangles
// E contains the boundary edges
void get_mesh(int n, Matrix &P, IntMatrix &T, IntMatrix &E)
{
	IntMatrix edges(4, 10000);
	int num_edges = 0;

	Matrix x(n + 1, 1);
	Matrix y(n + 1, 1);

	for(int i = 1; i <= n + 1; ++i)
	{
		x(i) = -1 + (i - 1) * 2. / n;
		y(i) = x(i);
	}

	P.resize(2, (n + 1) * (n + 1));
	for(int j = 1, count = 1; j <= n + 1; ++j)
	{
		for(int i = 1; i <= n + 1; ++i)
		{
			P(1, count) = x(i);
			P(2, count) = y(j);
			++count;
		}
	}

 	T.resize(10, 2 * n * n);
	for(int j = 1, count = 1; j <= n; ++j)
	{
		for(int i = 1; i <= n; ++i)
		{
			int index1 = (j - 1) * (n + 1) + i;
			int index2 = index1 + 1;
			int index3 = index2 + n + 1;
			int index4 = index3 - 1;
			T(1, count) = index1; T(2, count) = index2; T(3, count) = index3;
			++count;
			T(1, count) = index1; T(2, count) = index3; T(3, count) = index4;
			++count;
		}
	}

	int dof = 1;
	for(int i = 1; i <= T.get_num_cols(); ++i)
	{
		T(4, i) = dof;
		++dof;
	}

	for(int i = 1; i <= T.get_num_cols(); ++i)
	{
		//                        face1
		int index = find(num_edges, edges, T(2, i), T(3, i));
		if(index > 0)
		{
			T(5, i) = edges(4, index);
			T(6, i) = edges(3, index);
		}
		else
		{
			++num_edges;
			edges(1, num_edges) = T(2, i);
			edges(2, num_edges) = T(3, i);
			edges(3, num_edges) = T(5, i) = dof;
			edges(4, num_edges) = T(6, i) = dof + 1;
			dof += 2;
		}
		//                face2
		index = find(num_edges, edges, T(3, i), T(1, i));
		if(index > 0)
		{	
			T(7, i) = edges(4, index);
			T(8, i) = edges(3, index);
		}
		else
		{
			++num_edges;
			edges(1, num_edges) = T(3, i);
			edges(2, num_edges) = T(1, i);
			edges(3, num_edges) = T(7, i) = dof;
			edges(4, num_edges) = T(8, i) = dof + 1;
			dof += 2;
		}


		index = find(num_edges, edges, T(1, i), T(2, i));
		if(index > 0)
		{
			T(9, i) = edges(4, index);
			T(10, i) = edges(3, index);
		}
		else
		{
			++num_edges;
			edges(1, num_edges) = T(1, i);
			edges(2, num_edges) = T(2, i);
			edges(3, num_edges) = T(9, i) = dof;
			edges(4, num_edges) = T(10, i) = dof + 1;
			dof += 2;
		}
	}

	// configure E
	E.resize(2, 4 * n);
	for(int i = 1; i <= n; ++i)
	{
		E(1, i) = i;
		E(2, i) = i + 1;
	}

	for(int i = 1; i <= n; ++i)
	{
		E(1, n + i) = i * (n + 1);
		E(2, n + i) = (i + 1) * (n + 1);
	}

	for(int i = 1; i <= n; ++i)
	{
		E(1, 2 * n + i) = 1 + (i - 1) * (n + 1);
		E(2, 2 * n + i) = 1 + i * (n + 1);
	}

	for(int i = 1; i <= n; ++i)
	{
		E(1, 3 * n + i) = n * (n + 1) + i;
		E(2, 3 * n + i) = n * (n + 1) + i + 1;
	}
}
}//namespace easyfem


