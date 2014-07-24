/*
 * fem.h
 *
 *  Created on: Mar 3, 2014
 *      Author: lbj-lxy
 */

#ifndef FEM_H_
#define FEM_H_

#include <list>
#include "polynomial.h"

namespace easyfem
{
class Func
{
private:
	double coef[3][3];
public:
	Func()
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				coef[i][j] = 0.;
			}
		}
	}

	void add_monomial(int index1, int index2, double coef)
	{
		coef[index1][index2] += coef;
	}

	double eval(double x, double y)
	{
		double val = 0.;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				double tmp = coef[i][j];
				for (int k = 1; k <= i; ++k)
				{
					tmp *= x;
				}

				for (int k = 1; k <= j; ++k)
				{
					tmp *= y;
				}

				val += tmp;
			}
		}
		return val;
	}
};

class PKTriangle
{
public:
	static TMatrix<double> get_integration_info(int degree)
	{
		switch (degree)
		{
		case 1:
			TMatrix<double> info(4, 1);
			info(1, 1) = 1. / 3.;
			info(2, 1) = 1. / 3.;
			info(3, 1) = 1. / 3.;
			info(4, 1) = 1.;
			return info;
		case 2:
			TMatrix<double> info(4, 3);
			info(1, 1) = 2. / 3.;
			info(2, 1) = 1. / 6.;
			info(3, 1) = 1. / 6.;
			info(4, 1) = 1. / 3.;
			info(1, 2) = 1. / 6.;
			info(2, 2) = 2. / 3.;
			info(3, 2) = 1. / 6.;
			info(4, 2) = 1. / 3.;
			info(1, 3) = 1. / 6.;
			info(2, 3) = 1. / 6.;
			info(3, 3) = 2. / 3.;
			info(4, 3) = 1. / 3.;
			return info;
		case 3:
			TMatrix<double> info(4, 4);
			info(1, 1) = 1. / 3.;
			info(2, 1) = 1. / 3.;
			info(3, 1) = 1. / 3.;
			info(4, 1) = -27. / 48.;
			info(1, 2) = 0.6;
			info(2, 2) = 0.2;
			info(3, 2) = 0.2;
			info(4, 2) = 25. / 48.;
			info(1, 3) = 0.2;
			info(2, 3) = 0.6;
			info(3, 3) = 0.2;
			info(4, 3) = 25. / 48.;
			info(1, 4) = 0.2;
			info(2, 4) = 0.2;
			info(3, 4) = 0.6;
			info(4, 4) = 25. / 48.;
			return info;
		default:
			return TMatrix<double>(0, 0);
		}
	}
};

//              2
//             |  \
//             |   \
//             |    \
//             ------
//             3      1
class TriangleP1
{
public:
	static Func get_partial_x(int i)
	{
		Func func;
		if (i == 1)
		{
			func.add_monomial(0, 0, 1.);
		}
		else if (i == 2)
		{
			func.add_monomial(0, 0, 0.);
		}
		else
		{
			func.add_monomial(0, 0, -1.);
		}

		return func;
	}

	static Func get_partial_y(int i)
	{
		Func func;
		if (i == 1)
		{
			func.add_monomial(0, 0, 0.);
		}
		else if (i == 2)
		{
			func.add_monomial(0, 0, 1.);
		}
		else
		{
			func.add_monomial(0, 0, -1.);
		}

		return func;
	}
};
//                  2
//                 |  \
//                 6   5
//                 |    \
//                 3--4--1
class TriangleP2
{
public:
	static Func get_partial_x(int i)
	{
		Func func;
		if (i == 1)
		{
			func.add_monomial(0, 0, -1.);
			func.add_monomial(1, 0, 4.);
		}
		else if (i == 2)
		{
			func.add_monomial(0, 0, 0.);
		}
		else if (i == 3)
		{
			func.add_monomial(0, 0, -3.);
			func.add_monomial(1, 0, 4.);
			func.add_monomial(0, 1, 4.);
		}
		else if (i == 4)
		{
			func.add_monomial(0, 0, 4.);
			func.add_monomial(1, 0, -8.);
			func.add_monomial(0, 1, -4.);
		}
		else if (i == 5)
		{
			func.add_monomial(0, 1, 4.);
		}
		else if (i == 6)
		{
			func.add_monomial(0, 1, -4.);
		}

		return func;
	}

	static Func get_partial_y(int i)
	{
		Func func;
		if (i == 1)
		{
			func.add_monomial(0, 0, 0.);
		}
		else if (i == 2)
		{
			func.add_monomial(0, 0, -1.);
			func.add_monomial(0, 1, 4.);
		}
		else if (i == 3)
		{
			func.add_monomial(0, 0, -3.);
			func.add_monomial(1, 0, 4.);
			func.add_monomial(0, 1, 4.);
		}
		else if (i == 4)
		{
			func.add_monomial(1, 0, -4.);
		}
		else if (i == 5)
		{
			func.add_monomial(1, 0, 4.);
		}
		else if (i == 6)
		{
			func.add_monomial(0, 0, 4.);
			func.add_monomial(1, 0, -4.);
			func.add_monomial(0, 1, -8.);
		}

		return func;
	}

	static Func get_partial_xx(int i)
	{
		Func func;
		if (i == 1)
		{
			func.add_monomial(0, 0, 4.);
		}
		if (i == 3)
		{
			func.add_monomial(0, 0, 4.);
		}
		else if (i == 4)
		{
			func.add_monomial(0, 0, -8.);
		}

		return func;
	}

	static Func get_partial_xy(int i)
	{
		Func func;
		if (i == 3)
		{
			func.add_monomial(0, 0, 4.);
		}
		else if (i == 4)
		{
			func.add_monomial(0, 0, -4.);
		}
		else if (i == 5)
		{
			func.add_monomial(0, 0, 4.);
		}
		else if (i == 6)
		{
			func.add_monomial(0, 0, -4.);
		}

		return func;
	}

	static Func get_partial_yy(int i)
	{
		Func func;
		if (i == 2)
		{
			func.add_monomial(0, 0, 4.);
		}
		else if (i == 3)
		{
			func.add_monomial(0, 0, 4.);
		}
		else if (i == 6)
		{
			func.add_monomial(0, 0, -8.);
		}

		return func;
	}
};

//                  2
//                 |  \
//                 8    7
//                 |     \
//                 9      6
//                 |       \
//                 3--4--5--1
class TriangleP3
{
public:
	static Func get_partial_x(int i)
	{
		Func func;
		return func;
	}

	static Func get_partial_y(int i)
	{
		Func func;
		return func;
	}

	static Func get_partial_xx(int i)
	{
		Func func;
		return func;
	}

	static Func get_partial_xy(int i)
	{
		Func func;
		return func;
	}

	static Func get_partial_yy(int i)
	{
		Func func;
		return func;
	}
};
}//namespace easyfem

#endif /* FEM_H_ */
