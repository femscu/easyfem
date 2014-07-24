#include <cstdio>
#include "../INCLUDE/matrix.h"

namespace easyfem
{
//checked
static int find(int N, IntMatrix &S, int v1, int v2)
{
	for(int i = 1; i <= N; ++i)
	{
		if((S(1, i) == v1) && (S(2, i) == v2))
		{
			return i;
		}
		else if((S(1, i) == v2) && (S(2, i) == v1))
		{
			return i;
		}
	}

	return 0;
}

void number(IntMatrix &T)
{
	IntMatrix faces_found(4, 100000);
	int num_faces_found = 0;

	int num_elems = T.get_num_cols();
	int dof = 1;
	for(int i = 1; i <= num_elems; ++i)
	{
		T(4, i) = dof;
		++dof;
	}

	printf("num_elems = %d\n", num_elems);
	for(int i = 1; i <= num_elems; ++i)
	{
		//                             face1
		int index = find(num_faces_found, faces_found, T(2, i), T(3, i));
		if(index > 0)
		{
		 	T(5, i) = faces_found(4, index);
		 	T(6, i) = faces_found(3, index);
		}
		else
		{
			++num_faces_found;
			faces_found(1, num_faces_found) = T(2, i);
			faces_found(2, num_faces_found) = T(3, i);
			faces_found(3, num_faces_found) = T(5, i) = dof;
			faces_found(4, num_faces_found) = T(6, i) = dof + 1;
			dof += 2;
		}
		//                           face2
		index = find(num_faces_found, faces_found, T(3, i), T(1, i));
		if(index > 0)
		{
		 	T(7, i) = faces_found(4, index);
		 	T(8, i) = faces_found(3, index);
		}
		else
		{
			++num_faces_found;
			faces_found(1, num_faces_found) = T(3, i);
			faces_found(2, num_faces_found) = T(1, i);
			faces_found(3, num_faces_found) = T(7, i) = dof;
			faces_found(4, num_faces_found) = T(8, i) = dof + 1;
			dof += 2;
		}

		index = find(num_faces_found, faces_found, T(1, i), T(2, i));
		if(index > 0)
		{
		 	T(9, i) = faces_found(4, index);
		 	T(10, i) = faces_found(3, index);
		}
		else
		{
			++num_faces_found;
			faces_found(1, num_faces_found) = T(1, i);
			faces_found(2, num_faces_found) = T(2, i);
			faces_found(3, num_faces_found) = T(9, i) = dof;
			faces_found(4, num_faces_found) = T(10, i) = dof + 1;
			dof += 2;
		}
	}
}
void read_mesh(Matrix &P, IntMatrix &T)
{
	//*************************************************************
	//     read badic information of mesh
	//*************************************************************
	FILE *info_file = fopen("../../meshinfo.txt", "rt");
	if(info_file == 0)
	{
		printf("faile to open meshinfo.txt");
	}
	
	int num_nodes = 0;
	int num_elems = 0;
	double input = 0;
	fscanf(info_file, "%lf", &input);
	num_nodes = (int)input;
	fscanf(info_file, "%lf", &input);
	num_elems = (int)input;
	printf("num_nodes = %d\n", num_nodes);
	printf("num_elems = %d\n", num_elems);


	//****************************************************************
	//    read coordinates of nodes
	//***************************************************************
	FILE *node_file = fopen("../../node.txt", "rt");

	if (node_file == NULL)
	{
		printf("fail to open the file: node.txt");
		return;
	}

	P.resize(2, num_nodes);
	for(int i = 1; i <= num_nodes; ++i)
	{
		fscanf(node_file, "%lf", &P(1, i));
		fscanf(node_file, "%lf", &P(2, i));
	}


	//****************************************************************
	//     read information of elements
	//****************************************************************
	T.resize(10, num_elems);
	FILE *elem_file = fopen("../../elem.txt", "rt");
	if(elem_file == 0)
	{
		printf("faile to open elem.txt");
		return;
	}

	T.resize(10, num_elems);
	for(int i = 1; i <= num_elems; ++i)
	{
	 	double input;
		fscanf(elem_file, "%lf", &input);
		T(1, i) = input;
		fscanf(elem_file, "%lf", &input);
		T(2, i) = input;
		fscanf(elem_file, "%lf", &input);
		T(3, i) = input;
	}

	number(T);

	fclose(info_file);
	fclose(node_file);
	fclose(elem_file);
}

}// easyfem