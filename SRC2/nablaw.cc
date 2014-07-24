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
// compute the coeficents of \nabla_w{ the basis of W(\hat T)}
// \hat T is a triangle withe vetexes (0,0), (1, 0) and (0, 1).
//               3
//               | \
//               |  \
//               1---2
void nabla_w(Matrix &xx)
{
	xx.resize(6, 7);

	// mass matrix for P_1(\hat T)^2
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
	int num_triangles = triangles.get_num_cols();
	int dof = 0;
	for(int j = 1; j <= triangles.get_num_cols(); ++j)
	{
		for(int i = 4; i <= 10; ++i)
		{
			if(dof < triangles(i, j))
			{
				dof = triangles(i, j);
			}
		}
	}
	stiff_matrix.resize(dof, dof);

	// mass matrix for P_1(\hat T)^2
	Matrix mass(6, 6);
	mass(1, 1) = 1. / 12.; mass(1, 2) = 1. / 24.; mass(1, 3) = 1. / 24.; 
	mass(2, 1) = 1. / 24.; mass(2, 2) = 1. / 12.; mass(2, 3) = 1. / 24.;
	mass(3, 1) = 1. / 24.; mass(3, 2) = 1. / 24.; mass(3, 3) = 1. / 12.;
	mass(4, 4) = 1. / 12.; mass(4, 5) = 1. / 24.; mass(4, 6) = 1. / 24.;
	mass(5, 4) = 1. / 24.; mass(5, 5) = 1. / 12.; mass(5, 6) = 1. / 24.;
	mass(6, 4) = 1. / 24.; mass(6, 5) = 1. / 24.; mass(6, 6) = 1. / 12.;
		
	Matrix hatxx(0, 0);
	nabla_w(hatxx);
	
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
		inv(1, 1) =  a22 / det; inv(1, 2) = -a21 / det;
		inv(2, 1) = -a12 / det; inv(2, 2) =  a11 / det;

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
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 7; ++jj)
			{
				printf("%10f", xx(ii, jj));
			}
			printf("\n");
		}
		puts("-------------------------------------------");
		Matrix qu(7, 1);
		qu(1) = .125; qu(2) = .225; qu(3) = .275;
		qu(4) = .208333333; qu(5) = -.0416666667; qu(6) = -1. / 30.;
		qu(7) = 1. / 6.;
		Matrix lbj(6, 1);
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 7; ++jj)
			{
				lbj(ii) += xx(ii, jj) * qu(jj);
			}
		}

		for(int ii = 1; ii <= 6; ++ii)
		{
			printf("%f\n", lbj(ii));
		}
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
			stiff_matrix(index(ii), index(ii)) += area * product(6, xx.get_ptr(1, ii), xx.get_ptr(1, ii), mass);
			for(int jj = ii + 1; jj <= 7; ++jj)
			{
				stiff_matrix(index(ii), index(jj)) += area * product(6, xx.get_ptr(1, ii), xx.get_ptr(1, jj), mass);
				stiff_matrix(index(jj), index(ii)) = stiff_matrix(index(ii), index(jj));
			}
		}
	}
}
}//namespace easyfem