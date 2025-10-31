#include <iostream>
#include <iomanip>
#include "Matrix.h"
using namespace std;
using namespace lstd;
int main() {
	cout << fixed << std::setprecision(2);
	int n;
	cin >> n;
	Matrix co_matrix(n, n);
	Matrix col_vec(n, 1);
	cin >> co_matrix >> col_vec;
	if (co_matrix.isSingular()) cout << "Error";
	else cout << (co_matrix ^ (-1)) * col_vec;
}