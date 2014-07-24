#include <vector>
#include <cstdio>
#include <cmath>

#include "../INCLUDE/matrix.h"

namespace easyfem
{

double norm(double x, double y)
{
	return sqrt(x * x + y * y);
}
double inner_product(Matrix &x, Matrix &y)
{
	double ret = 0.;

	int len = x.size();
	for(int i = 1; i <= len; ++i)
	{
		ret += x(i) * y(i);
	}

	return ret;
}
double inner_product(Matrix &x, Matrix &y, Matrix &A)
{
	int len = x.size();
	Matrix tmp(len, 1);
	for(int i = 1; i <= len; ++i)
	{
		for(int j = 1; j <= len; ++j)
		{
			tmp(i) += A(i, j) * x(j);
		}
	}

	return inner_product(tmp, y);
}
double inner_product(double *x, double *y, Matrix &A)
{
	int n = A.get_rows();
	Matrix tmp(n, 1);
	for(int i = 1; i <= n; ++i)
	{
		for(int j = 1; j <= n; ++j)
		{
			tmp(i) += A(i, j) * x[j - 1];
		}
	}

	double ret = 0.;
	for(int i = 1; i <= n; ++i)
	{
		ret += tmp(i) * y[i - 1];
	}

	return ret;
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
void calc_source_vec(Matrix &nodes, IntMatrix &triangles, Matrix &source)
{
	source.resize(triangles.get_cols(), 1);
	for(int i = 1; i <= triangles.get_cols(); ++i)
	{
		double x1 = nodes(1, triangles(1, i));
		double y1 = nodes(2, triangles(1, i));

		double x2 = nodes(1, triangles(2, i));
		double y2 = nodes(2, triangles(2, i));

		double x3 = nodes(1, triangles(3, i));
		double y3 = nodes(2, triangles(3, i));

		double area = 0.5 * (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
		source(i) = area;
	}
}
// k = 0
void calc_stiff_matrix(Matrix &nodes, IntMatrix &triangles, Matrix &stiff_matrix)
{
	int num_triangles = triangles.get_cols();

	// local mass matrix for P_1(T)^2
	Matrix local_mass(6, 6);
	local_mass.set_zeros();
	for(int i = 1; i <= 6; ++i)
	{
		local_mass(i, i) = 1. / 6.;
	}
	local_mass(1, 2) = 1. / 12.; local_mass(1, 3) = 1. / 12.;
	local_mass(2, 1) = 1. / 12.; local_mass(2, 3) = 1. / 12.;
	local_mass(3, 1) = 1. / 12.; local_mass(3, 2) = 1. / 12.;
	local_mass(4, 5) = 1. / 12.; local_mass(4, 6) = 1. / 12.;
	local_mass(5, 4) = 1. / 12.; local_mass(5, 6) = 1. / 12.;
	local_mass(6, 4) = 1. / 12.; local_mass(6, 5) = 1. / 12.;
	// end local mass matrix of P_1(T)^2

	// the inverse of local mass matrix of P_1(T)^2
	Matrix inv_local_mass(6, 6);
	inv_local_mass.set_zeros();
	for(int i = 1; i <= 6; ++i)
	{
		inv_local_mass(i, i) = 9.;
	}
	inv_local_mass(1, 2) = -3.; inv_local_mass(1, 3) = -3.;
	inv_local_mass(2, 1) = -3.; inv_local_mass(2, 3) = -3.;
	inv_local_mass(3, 1) = -3.; inv_local_mass(3, 2) = -3.;
	inv_local_mass(4, 5) = -3.; inv_local_mass(4, 6) = -3.;
	inv_local_mass(5, 4) = -3.; inv_local_mass(5, 6) = -3.;
	inv_local_mass(6, 4) = -3.; inv_local_mass(6, 5) = -3.;
	//end inverse of local mass matrix of P_1^(T)^2

	// degree of freedoms
	int dof = 7 * triangles.get_cols();

	// stiffness matrix
	stiff_matrix.resize(dof, dof);
	stiff_matrix.set_zeros();
	for(int i = 1; i <= num_triangles; ++i)
	{
		double x1 = nodes(1, triangles(1, i));
		double y1 = nodes(2, triangles(1, i));

		double x2 = nodes(1, triangles(2, i));
		double y2 = nodes(2, triangles(2, i));

		double x3 = nodes(1, triangles(3, i));
		double y3 = nodes(2, triangles(3, i));

		double area = 0.5 * (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
		double area_face1 = sqrt((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3));
		double area_face2 = sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3));
		double area_face3 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

		printf("area = %f\n", area);
		// nabla_lambda(1) is \partial_x\lambda_1, nabla_lambda(2) is \partial_y\lambda1
		// nabla_lambda(3) is \partial_x\lambda_2, nabla_lambda(4) is \partial_y\lambda2
		// ...
		Matrix nabla_lambda(6, 1);
		nabla_lambda(1) = y2 - y3; nabla_lambda(2) = x3 - x2;
		nabla_lambda(3) = y3 - y1; nabla_lambda(4) = x1 - x3;
		nabla_lambda(5) = y1 - y2; nabla_lambda(6) = x2 - x1;
		//*********************************************************************
		//
		//*********************************************************************
		Matrix xx(6, 7);
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 6; ++jj)
			{
				xx(ii, 1) -= .5 * inv_local_mass(ii, jj) * nabla_lambda(jj);
			}
			xx(ii, 1) /= area;
		}

		// the outward normal vectors of three faces
		Matrix normal_vectors(2, 3);
		normal_vectors(1, 1) = (y3 - y2) / area_face1;
		normal_vectors(2, 1) = (x2 - x3) / area_face1;
		normal_vectors(1, 2) = (y1 - y3) / area_face2;
		normal_vectors(2, 2) = (x3 - x1) / area_face2;
		normal_vectors(1, 3) = (y2 - y1) / area_face3;
		normal_vectors(2, 3) = (x1 - x2) / area_face3;

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//                            1--face1
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		Matrix b(6, 1);
		b.set_zeros();
		b(2) = area_face1 / 3. * normal_vectors(1, 1);
		b(3) = area_face1 / 6. * normal_vectors(1, 1);
		b(5) = area_face1 / 3. * normal_vectors(2, 1);
		b(6) = area_face1 / 6. * normal_vectors(2, 1);

		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 6; ++jj)
			{
				xx(ii, 2) += inv_local_mass(ii, jj) * b(jj);
			}
			xx(ii, 2) /= area;
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//                             2--face1
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		b.set_zeros();
		b(2) += area_face1 / 6. * normal_vectors(1, 1);
		b(3) += area_face1 / 3. * normal_vectors(1, 1);
		b(5) += area_face1 / 6. * normal_vectors(2, 1);
		b(6) += area_face1 / 3. * normal_vectors(2, 1);

		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 6; ++jj)
			{
				xx(ii, 3) += inv_local_mass(ii, jj) * b(jj);
			}
			xx(ii, 3) /= area;
		}
	    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//                              3--face2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		b.set_zeros();
		b(1) += area_face2 / 6. * normal_vectors(1, 2);
		b(3) += area_face2 / 3. * normal_vectors(1, 2);
		b(4) += area_face2 / 6. * normal_vectors(2, 2);
		b(6) += area_face2 / 3. * normal_vectors(2, 2);
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 6; ++jj)
			{
				xx(ii, 4) = inv_local_mass(ii, jj) * b(jj);
			}
			xx(ii, 4) /= area;
		}
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//                               4--face2
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		b.set_zeros();
		b(1) += area_face2 / 3. * normal_vectors(1, 2);
		b(3) += area_face2 / 6. * normal_vectors(1, 2);
		b(4) += area_face2 / 3. * normal_vectors(2, 2);
		b(6) += area_face2 / 6. * normal_vectors(2, 2);
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 6; ++jj)
			{
				xx(ii, 5) += inv_local_mass(ii, jj) * b(jj);
			}
			xx(ii, 5) /= area;
		}
		//**********************************************************************
		//                               5--face3
		//**********************************************************************
		b.set_zeros();
		b(1) += area_face3 / 3. * normal_vectors(1, 3);
		b(2) += area_face3 / 6. * normal_vectors(1, 3);
		b(4) += area_face3 / 3. * normal_vectors(2, 3);
		b(5) += area_face3 / 6. * normal_vectors(2, 3);
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 6; ++jj)
			{
				xx(ii, 6) += inv_local_mass(ii, jj) * b(jj);
			}
			xx(ii, 6) /= area;
		}
		//**********************************************************************
		//                               6--face3
		//**********************************************************************
		b.set_zeros();
		b(1) += area_face3 / 6. * normal_vectors(1, 3);
		b(2) += area_face3 / 3. * normal_vectors(1, 3);
		b(4) += area_face3 / 6. * normal_vectors(2, 3);
		b(5) += area_face3 / 3. * normal_vectors(2, 3);
		for(int ii = 1; ii <= 6; ++ii)
		{
			for(int jj = 1; jj <= 6; ++jj)
			{
				xx(ii, 7) += inv_local_mass(ii, jj) * b(jj);
			}
			xx(ii, 7) /= area;
		}

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
			stiff_matrix(index(ii), index(ii)) = area * inner_product(xx.get_ptr(1, ii), xx.get_ptr(1, ii), local_mass);
			for(int jj = ii + 1; jj <= 7; ++jj)
			{
				stiff_matrix(index(ii), index(jj)) += area * inner_product(xx.get_ptr(1, ii), xx.get_ptr(1, jj), local_mass);
				stiff_matrix(index(jj), index(ii)) = stiff_matrix(index(ii), index(jj));
			}
		}
	}
}

}//namespace easyfem

