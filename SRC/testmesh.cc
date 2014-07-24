#include <cstdio>
#include <vector>
#include "../INCLUDE/matrix.h"

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
	Matrix nodes(0, 0);
	IntMatrix triangles(0, 0);

	get_mesh(2, nodes, triangles);

	for(int i = 1; i <= triangles.get_cols(); ++i)
	{
		printf("%d:", i);
		for(int j = 1; j <= 10; ++j)
		{
			printf("%d,", triangles(j, i));
		}
		printf("\n");
	}
	SaveData("/home/lxy/workspace/nodes.text", nodes.get_cols() * 2, nodes.get_ptr(1));
	SaveData("/home/lxy/workspace/T.text", triangles.get_cols() * 10, triangles.get_ptr(1));
	printf("-----------------------------Done!----------------\n");
	return 0;
}