#include <cstdio>
#include <vector>
#include "../INCLUDE/matrix.h"

namespace easyfem
{

static double product(int N, double *X, double *Y, Matrix &A)
{
	Matrix Z(N, 1);
	for(int i = 1; i <= N; ++i)
	{
		double tmp = 0;
		for(int j = 1; j <= N; ++j)
		{
			tmp += A(i, j) * X[j - 1];
		}
		Z(i) = tmp;
	}

	double ret = 0.;
	for(int i = 1; i <= N; ++i)
	{
		ret += Z(i) * Y[i - 1];
	}

	return ret;
}
// T is a triangle withe vetexes (0,0), (1, 0) and (0, 1).
//               3
//               | \
//               |  \
//               1---2
void nabla_w(Matrix &xx)
{
	xx.resize(6, 7);

	// mass matrix for P_1(T)^2
	Matrix mass(6, 6);
	mass(1, 1) = 1. / 12.; mass(1, 2) = 1. / 24.; mass(1, 3) = 1. / 24.; 
	mass(2, 1) = 1. / 24.; mass(2, 2) = 1. / 12.; mass(2, 3) = 1. / 24.;
	mass(3, 1) = 1. / 24.; mass(3, 2) = 1. / 24.; mass(3, 3) = 1. / 12.;
	mass(4, 4) = 1. / 12.; mass(4, 5) = 1. / 24.; mass(4, 6) = 1. / 24.;
	mass(5, 4) = 1. / 24.; mass(5, 5) = 1. / 12.; mass(5, 6) = 1. / 24.;
	mass(6, 4) = 1. / 24.; mass(6, 5) = 1. / 24.; mass(6, 6) = 1. / 12.;
	
	// the inverse of mass
	Matrix inv_mass(6, 6);
	inv_mass(1, 1) = 18.; inv_mass(1, 2) = -6.; inv_mass(1, 3) = -6.; 
	inv_mass(2, 1) = -6.; inv_mass(2, 2) = 18.; inv_mass(2, 3) = -6.;
	inv_mass(3, 1) = -6.; inv_mass(3, 2) = -6.; inv_mass(3, 3) = 18.;
	inv_mass(4, 4) = 18.; inv_mass(4, 5) = -6.; inv_mass(4, 6) = -6.;
	inv_mass(5, 4) = -6.; inv_mass(5, 5) = 18.; inv_mass(5, 6) = -6.;
	inv_mass(6, 4) = -6.; inv_mass(6, 5) = -6.; inv_mass(6, 6) = 18.;

	Matrix b(6, 7);
	//                   1
	b(1, 1) = .5; b(2, 1) = -.5; b(4, 1) = .5; ; b(6, 1) = -.5;
	//                   2                               
	b(2, 2) = b(5, 2) = 1. / 3.;
	b(3, 2) = b(6, 2) = 1. / 6.;

	//                   3
	b(2, 3) = b(5, 3) = 1. / 6.;
	b(3, 3) = b(6, 3) = 1. / 3.;

	//                   4
	b(1, 4) = -1. / 6.;
	b(3, 4) =  -1. / 3.;
	//                   5
	b(1, 5) = -1. / 3.;
	b(3, 5) = -1. / 6.;
	//                   6
	b(4, 6) = -1. / 3.;
	b(5, 6) = -1. / 6.;
	//                   7
	b(4, 7) = -1. / 6.;
	b(5, 7) = -1. / 3.;

	for(int j = 1; j <= 7; ++j)
	{
		for(int i = 1; i <= 6; ++i)
		{
			double tmp = 0;
			for(int k = 1; k <= 6; ++k)
			{
				tmp += inv_mass(i, k) * b(k, j);
			}
			xx(i, j) = tmp;
		}
	}

	return;
}



void calc_stiff_matrix(Matrix &nodes, IntMatrix &triangles, Matrix &stiff_matrix)
{
	int num_triangles = triangles.get_cols();
	int dof = 7 * num_triangles;
	stiff_matrix.resize(dof, dof);

	// mass matrix for P_1(T)^2
	Matrix mass(6, 6);
	mass(1, 1) = 1. / 12.; mass(1, 2) = 1. / 24.; mass(1, 3) = 1. / 24.; 
	mass(2, 1) = 1. / 24.; mass(2, 2) = 1. / 12.; mass(2, 3) = 1. / 24.;
	mass(3, 1) = 1. / 24.; mass(3, 2) = 1. / 24.; mass(3, 3) = 1. / 12.;
	mass(4, 4) = 1. / 12.; mass(4, 5) = 1. / 24.; mass(4, 6) = 1. / 24.;
	mass(5, 4) = 1. / 24.; mass(5, 5) = 1. / 12.; mass(5, 6) = 1. / 24.;
	mass(6, 4) = 1. / 24.; mass(6, 5) = 1. / 24.; mass(6, 6) = 1. / 12.;
		
	Matrix hatxx(0, 0);
	nabla_w(hatxx);
	//*
		puts("--------------hatxx--------------");
		for(int i = 1; i <= 6; ++i)
		{
			for(int j = 1; j <= 7; ++j)
			{
				printf("%f,", hatxx(i, j));
			}
			printf("\n");
		}
		puts("------------end hatxx----------------");
	//*/
	Matrix xx(6, 7);
	Matrix inv(2, 2);
	for(int i = 1; i <= num_triangles; ++i)
	{
		double x1 = nodes(1, triangles(1, i));
		double y1 = nodes(2, triangles(1, i));
		double x2 = nodes(1, triangles(2, i));
		double y2 = nodes(2, triangles(2, i));
		double x3 = nodes(1, triangles(3, i));
		double y3 = nodes(2, triangles(3, i));

		//   (a11 a12)
		//   (a21 a22)   is maps \hat T to T
		double a11 = x2 - x1; 
		double a12 = x3 - x1;
		double a21 = y2 - y1;
		double a22 = y3 - y1;
		double det = a11 * a22 - a12 * a21;
		// inv := (a11 a12)^{-T}
		//        (a21 a22)
		inv(1, 1) = a22 / det; inv(1, 2) = -a21 / det;
		inv(2, 1) = -a12 / det; inv(2, 2) = a11 / det;

		//   \hat{\nabla_w v} 	= 	A^{-T} \hat\nabla_w\hat v
		// \hat{\nabla_w \mu} 	= 	A^{-T} \hat\nabla_w\hat\mu
		xx.set_zeros();
		for(int jj = 1; jj <= 7; ++jj)
		{
			for(int ii = 1; ii <= 3; ++ii)
			{
				xx(ii, jj) += inv(1, 1) * hatxx(ii, jj);
				xx(ii, jj) += inv(1, 2) * hatxx(ii + 3, jj);
			}

			for(int ii = 4; ii <= 6; ++ii)
			{
				xx(ii, jj) += inv(2, 1) * hatxx(ii - 3, jj);
				xx(ii, jj) += inv(2, 2) * hatxx(ii, jj); 
			}
		}
		//*
		puts("--------------------xx---------------");
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 7; ++jj)
			{
				printf("%f,", xx(ii, jj));
			}
			printf("\n");
		}
		puts("-----------------end xx-----------------");
		//*/
		double area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
		//**********************************************************************
		//                       Assemble local stiffness matrix
		//**********************************************************************
		IntMatrix index(7,1);
		for(int ii = 1; ii <= 7; ++ii)
		{
			index(ii) = triangles(3 + ii, i);
		}

		for(int ii = 1; ii <= 7; ++ii)
		{
			stiff_matrix(index(ii), index(ii)) = area * product(6, xx.get_ptr(1, ii), xx.get_ptr(1, ii), mass);
			for(int jj = ii + 1; jj <= 7; ++jj)
			{
				stiff_matrix(index(ii), index(jj)) += area * product(6, xx.get_ptr(1, ii), xx.get_ptr(1, jj), mass);
				stiff_matrix(index(jj), index(ii)) = stiff_matrix(index(ii), index(jj));
			}
		}
	}
}


bool is_bnd(Matrix &nodes, int index)
{
	if(nodes(1, index) == -1 || nodes(1, index) == 1 || nodes(2, index) == -1 || nodes(2, index) == 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool check2(Matrix &nodes, int index1, int index2)
{
	double x1 = nodes(1, index1); double y1 = nodes(2, index1);
	double x2 = nodes(1, index2); double y2 = nodes(2, index2);
	if(x1 == x2)
	{
		return true;
	}
	if(y1 == y2)
	{
		return true;
	}

	return false;
}
void calc_bound_con(Matrix &nodes, IntMatrix &triangles, std::vector<int> &bnd)
{
	bnd.resize(1000);
	int count = 0;
	for(int i = 1; i <= triangles.get_cols(); ++i)
	{
		printf("count = %d\n", count);
		bool vert1 = is_bnd(nodes, triangles(1, i));
		bool vert2 = is_bnd(nodes, triangles(2, i));
		bool vert3 = is_bnd(nodes, triangles(3, i));

		printf("%d, %d, %d\n", vert1, vert2, vert3);
		printf("%d\n", triangles.get_rows());
		if(vert1 && vert2)
		{
			if(check2(nodes, triangles(1, i), triangles(2, i)))
			{
				bnd[count] = triangles(9, i);
				bnd[count + 1] = triangles(10, i);
				count += 2;
			}
		}

		if(vert1 && vert3)
		{
			if(check2(nodes, triangles(1, i), triangles(3, i)))
			{
				bnd[count] = triangles(7, i);
				bnd[count + 1] = triangles(8, i);
				count += 2;
			}
		}

		if(vert2 && vert3)
		{
			if(check2(nodes, triangles(2, i), triangles(3, i)))
			{
				bnd[count] = triangles(5, i);
				bnd[count + 1] = triangles(6, i);
				count += 2;
			}
		}
	}

	bnd.resize(count);
}

}//namespace easyfem