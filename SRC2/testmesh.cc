#include <cstdio>
#include "../INCLUDE/matrix.h"

namespace easyfem
{
void get_mesh(int n, Matrix &P, IntMatrix &T, IntMatrix &E);
}

using namespace easyfem;

int main()
{
	Matrix P(0, 0);
	IntMatrix T(0, 0);
	IntMatrix E(0, 0);

	get_mesh(3, P, T, E);

	for(int i = 1; i <= T.get_num_cols(); ++i)
	{
		printf("%d:", i);
		for(int j = 1; j <= 10; ++j)
		{
			printf("%d,", T(j, i));
		}
		printf("\n");
	}

	for(int j = 1; j <= E.get_num_cols(); ++j)
	{
		printf("%d----%d\n", E(1, j), E(2, j));
	}

	printf("-----------------------------Done!----------------\n");
	return 0;
}