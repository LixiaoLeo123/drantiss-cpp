/*------------------------
* |      cpp utils       |
* |     by Drantiss      |
* |         2025         |
* ------------------------
* IMPORTANT: All matrix indices are 1-based!
*/
#include <iostream>
#include <iomanip>
#include "Matrix.h"
using namespace std;
using namespace drantiss;
int main() {
	//Solving linear equations with Matrix.h
	cout << fixed << std::setprecision(2);
	int n;
	cin >> n;
	Matrix co_matrix(n, n);
	Matrix col_vec(n, 1);
	cin >> co_matrix >> col_vec;
	if (co_matrix.isSingular()) cout << "Error";
	else cout << (co_matrix ^ (-1)) * col_vec;
	//Basic matrix operations
	Matrix addition = co_matrix + co_matrix ^ (-2);
	addition *= ~co_matrix *= Matrix::I(n);
	addition.addColMultiple(1, 4, 2);
	if (n > 4) {
		addition = addition.getMinor(3, 4);
		double cofactor = addition.cofactor(3, 4);
	}
	//Some methods for constructing and initializing matrices
	Matrix matrixWithI = Matrix::I(n);
	Matrix matrixWithDiag = Matrix::diag(new double[5] {1, 1, 2, 3, 5}, 5);   //No use of mobile semantics
	Matrix vectorWithColvec = Matrix::colvec(new double[4] {1, 2, 3, 4}, 4);  //Attention: It may result in reading out of heap memory
	Matrix newMatrix(new double[6] {1, 1, 2, 2, 3, 5}, 2, 3);
	//See Matrix.h for more usage
}