/*
 * matrix.h
 *
 *  Created on: Mar 18, 2014
 *      Author: libinjie
 */

#ifndef MATRIX_H_
#define MATRIX_H_
#include <cstdio>
namespace easyfem
{
template<typename T>
class TMatrix
{
private:
	int num_rows;
	int num_cols;
	T *ptr;
	//auxiliary pointers
	T *ptr_matrix;
	T *ptr_vector;
public:
	struct Slice
	{ 
		int num_rows;
		int num_cols;
		int num_rows_parent;
		T *ptr;
		// auxiliary pointers
		T *ptr_matrix;
		
		Slice(TMatrix& parent, int index1, int index2, int nrows, int ncols)
		{
			num_rows = nrows;
			num_cols = ncols;
			num_rows_parent = parent.get_num_rows();
			ptr = parent.get_ptr(index1, index2);
			ptr_matrix = ptr -1 - num_row_parent;
		}

		T &operator()(int i, int j)
		{
			return *(ptr_matrix + i + j * num_rows_parent); 	
		}
	};	
public:
	TMatrix(int nrows, int ncols):
		num_rows(nrows), num_cols(ncols),
		ptr(0), ptr_matrix(0), ptr_vector(0)
	{
		ptr = new T[nrows * ncols];
		ptr_matrix = ptr - 1 -num_rows;
		ptr_vector = ptr - 1;

		this->set_zeros();
	}

	int get_num_rows() const
	{
		return num_rows;
	}

	int get_num_cols() const
	{
		return num_cols;
	}

	T *get_ptr()
	{
		return ptr;
	}

	T *get_ptr(int i, int j)
	{
		if(i < 1 || i > num_rows || j < 1 || j > num_cols)
			throw "out of range";

		return ptr_matrix + i + j * num_rows;
	}

	T *get_ptr(int i)
	{
		if(i < 1 || i > num_rows * num_cols)
			throw "out of range";

		return ptr_vector + i;
	}

	void resize(int nrows, int ncols)
	{

		if(ptr)
		{
			delete []ptr;
		}
		
		ptr = new T[nrows * ncols];
		num_rows = nrows;
		num_cols = ncols;
		ptr_matrix = ptr - 1 - num_rows;
		ptr_vector = ptr - 1;

		this->set_zeros();
	}

	int size() const
	{
		return num_rows * num_cols;
	}
	void swap(TMatrix &lhs, TMatrix &rhs)
	{
		int tmp = rhs.num_rows;
		rhs.num_rows = lhs.num_rows;
		lhs.num_rows = tmp;

		tmp = rhs.num_cols;
		rhs.num_cols = lhs.num_cols;
		lhs.num_cols = tmp;

		//******************************
		T *tmp_ptr = rhs.ptr;
		rhs.ptr = lhs.ptr;
		lhs.ptr = tmp_ptr;

		tmp_ptr = rhs.ptr_matrix;
		rhs.ptr_matrix = lhs.ptr_matrix;
		lhs.ptr_matrix = tmp_ptr;

		tmp_ptr = rhs.ptr_vector;
		rhs.ptr_vector = lhs.ptr_vector;
		lhs.ptr_vector = tmp_ptr;
	}
	void set_zeros()
	{
		for(int i = 0; i < num_rows * num_cols; ++i)
		{
			ptr[i] = T(0);
		}
	}

	T &operator()(int i,int j)
	{
		if(i < 1 || i > num_rows || j < 1 || j > num_cols)
			throw "out of range";

		return ptr_matrix[i + j * num_rows];
	}

	T &operator()(int i)
	{
		if(i < 1 || i > num_rows * num_cols)
			throw "out of range";

		return ptr_vector[i];
	}
}; 

typedef TMatrix<int> IntMatrix;
typedef TMatrix<double> Matrix;
}//namespace easyfem 

#endif /* MATRIX_H_ */
