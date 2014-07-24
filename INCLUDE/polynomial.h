///*
// * polynomial.h
// *
// *  Created on: Feb 28, 2014
// *      Author: lbj-lxy
// */
//
//#ifndef POLYNOMIAL_H_
//#define POLYNOMIAL_H_
//
//namespace easyfem
//{
//int  factorial(int n)
//{
//	if(n == 1)
//	{
//		return 1;
//	}
//	else
//	{
//		return n * factorial(n - 1);
//	}
//}
//
//int nchoosek(int n, int m)
//{
//	n = n < m ? m : n;
//
//	if(n == m)
//	{
//		return 1;
//	}
//	else
//	{
//		return factorial(n) / factorial(m) / factorial(n - m);
//	}
//}
//
//class MultiIndex
//{
//private:
//	int dim;
//	IntMatrix mindex;
//public:
//	MultiIndex(int dim)
//	{
//		this->dim = dim;
//		mindex.resize(dim, 1);
//		mindex.set_zeros();
//	}
//
//	MultiIndex(int dim, IntMatrix &mindex)
//	{
//		this->dim = dim;
//		this->mindex = mindex;
//	}
//
//	int get_degree()
//	{
//		return MultiIndex::get_degree(mindex.get_conslice(2, mindex.get_size()));
//	}
//
//	static int get_degree(const IntMatrix::ConSlice &mindex)
//	{
//		int degree = 0;
//		for(int i = 1; i <= mindex.get_size(); ++i)
//		{
//			degree += mindex(i);
//		}
//
//		return degree;
//	}
//
//	static int calc_global_index(int dim, const IntMatrix::ConSlice &mindex)
//	{
//		if(dim == 1)
//		{
//			return mindex(1) + 1;
//		}
//
//		int degree = MultiIndex::get_degree(mindex);
//		if(degree == 0)
//		{
//			return 1;
//		}
//
//		int index = 0;
//		index += nchoosek(dim + degree - 1, dim);
//		return index + MultiIndex::calc_global_index(dim - 1, mindex.get_conslice(2, mindex.get_size()));
//	}
//
//	int get_global_index() const
//	{
//		return MultiIndex::calc_global_index(dim, mindex.get_conslice(1, mindex.get_size()));
//	}
//
////	int operator()(int i) const
//	{
//		return mindex(i);
//	}
//};
//
//
//template<typename T>
//class Polynomial
//{
//private:
//	int dim;
//	DynamicArray<T, 32> coefs;
//public:
//	Polynomial();
//
//	int get_degree() const
//	{
//		int index_last_nonzero = coefs.get_size();
//		while((index_last_nonzero > 0) && (coefs(index_last_nonzero) == T(0)))
//		{
//			--index_last_nonzero;
//		}
//
//		if(index_last_nonzero == 0)
//		{
//			return -1;
//		}
//		else
//		{
//
//		}
//
//	}
//
//	void add_monomial(const T &coef, const MultiIndex &mindex)
//	{
//		int index = mindex.get_global_index();
//		coefs(index) += coef;
//	}
//
//
//};
//} // namespace
//
//
//
//
//#endif /* POLYNOMIAL_H_ */
