My own implementation of Compressed sparse row format for sparse matrix storage: https://en.wikipedia.org/wiki/Sparse_matrix


## Clone 
```
git clone https://github.com/all4one1/FromOuterSparse.git
```

## Example usage 

```cpp
#include "FromOuterSparse/SparseMatrix.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
using std::cout;
using std::endl;


using Matrix = std::vector<std::vector<double>>;

int main()
{
	Matrix example =
	{ { 30,3,4,0,0,0 },
	 { 4,22,1,3,0,0 },
	 { 5,7,33,6,7,0 },
	 { 0,1,2,42,3,3 },
	 { 0,0,2,11,52,2 },
	 { 0,0,0,3,9,26 } };

	SparseMatrix SM;
	SM = example;

	SM.print_compressed_matrix();

	//output:
	//30 3 4 4 22 1 3 5 7 33 6 7 1 2 42 3 3 2 11 52 2 3 9 26
	//0 1 2 0 1 2 3 0 1 2 3 4 1 2 3 4 5 2 3 4 5 3 4 5
	//0 3 7 12 17 21 24

	//access an element SM[1][1] = 22
	cout << SM[1][1] << endl; // lvalue, produces 0 if doesn't exist
	cout << SM(1, 1) << endl; // lvalue, produces 0 if doesn't exist
	cout << SM.get_element(1, 1) << endl; // rvalue, save access


	//save or recover a sparse matrix
	SM.save_compressed_matrix();
	SM.read_compressed_matrix();


	return 0;
}
```
