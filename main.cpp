#include "SparseMatrix.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

using std::cout;
using std::endl;

int main()
{
	//Let us make a test
	SparseMatrix SM(6);
	double A[6][6] =
	{
		{ 30,3,4,0,0,0 },
		{ 4,22,1,3,0,0 },
		{ 5,7,33,6,7,0 },
		{ 0,1,2,42,3,3 },
		{ 0,0,2,11,52,2 },
		{ 0,0,0,3,9,26 },
	};

	double** M = new double* [6];
	for (int i = 0; i < 6; i++)
	{
		M[i] = new double[6]; 
	}

	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			M[i][j] = A[i][j]; 
		}
	}

	SM.make_sparse_from_double_array(6, M);
	SM.print_all();	
	SM.print_sequently();

	//access an element
	cout << SM[1][1] << endl; // lvalue, produce 0 if doesn't exists
	cout << SM(1, 1) << endl; // lvalue, produce 0 if doesn't exists
	cout << SM.get_element(1, 1) << endl; // rvalue, save access

	cout << "End" << endl;
	return 0;
}