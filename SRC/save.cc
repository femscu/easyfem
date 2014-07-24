/*
 * save.cpp
 *
 *  Created on: 2013-8-20
 *      Author: Li Binjie
 */

#include<cstdio>

namespace easyfem
{
typedef float real;
typedef double doublereal;
template<typename T>
void SaveData(const char *filename, const int N, const T *X)
{
	FILE *fp = fopen(filename, "wt");

	if (fp == NULL)
	{
		printf("fail to create the file: %s\n", filename);
		return;
	}

	for (int i = 0; i < N; ++i)
	{
		fprintf(fp, "%.16f\t", X[i]);
	}

	fclose(fp);
}

void SaveData(const char *filename, const int N, const int *X)
{
	FILE *fp = fopen(filename, "wt");

	if (fp == NULL)
	{
		printf("fail to create the file: %s\n", filename);
		return;
	}

	for (int i = 0; i < N; ++i)
	{
		fprintf(fp, "%d\t", X[i]);
	}

	fclose(fp);
}
void SaveData(const char *filename, const int N, const double *X)
{
	SaveData<double>(filename, N, X);
}


} // namespace scu
