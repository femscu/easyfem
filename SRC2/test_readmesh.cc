#include <cstdio>
#include "../INCLUDE/matrix.h"

namespace easyfem
{
void read_mesh(Matrix &P, IntMatrix &T);
}

using namespace easyfem;

int main()
{
	Matrix P(0, 0);
	IntMatrix T(0, 0);
	
	read_mesh(P, T);

	for(int i = 1; i <= P.get_num_cols(); ++i)
	{
		printf("%d: %f, %f, %f\n", i, P(1, i), P(2, i), P(3, i));
	}

	puts("------------------------------T-------------");
	for(int i = 1; i <= T.get_num_cols(); ++i)
	{
		printf("%d:", i);
		for(int j = 1; j <= 10; ++j)
		{
			printf("%3d", T(j, i));
		}
		printf("\n");
	}
	return 0;
}