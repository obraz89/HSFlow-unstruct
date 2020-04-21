#include "matrix_small.h"

#include "logging.h"

//#include <cmath>

t_Vec3 operator*(double val, const t_Vec3& vec) {
	t_Vec3 ret = vec;
	ret *= val;
	return ret;
};

// ******************************t_MatRotN

// cos_psi*cos_tet = nx
// sin_psi*cos_tet = ny
// sin_tet = -nz
void t_MatRotN::calc_rot_angles_by_N(const t_Vec3& N) {

	const double& nx = N[0];
	const double& ny = N[1];
	const double& nz = N[2];

	// if Nz~1, psi is arbitrary, choose psi=0
	// also use empirical tolerance to detect this situation
	// if epsilon is machine eps for a given arithmetics: 1+epsilon = 1
	// (for double precision epsilon is 10^-16), then minimum tolerance is sqrt(eps)~10^-8
	double tol = 1.0e-07;
	if (fabs(fabs(nz) - 1.0) < tol) {

		cos_psi = 1.0;
		sin_psi = 0.0;

		sin_tet = -1.0*sign(nz);
		cos_tet = 0.0;


	}
	else {

		sin_tet = -1.0 * nz;
		// convention: choose cos_tet>0
		cos_tet = sqrt(1.0-sin_tet*sin_tet);

		cos_psi = nx / cos_tet;
		sin_psi = ny / cos_tet;

	}

}
// ******************************t_SqMat3
// direct rotation matrix given by angles 
void t_SqMat3::set(const t_MatRotN& rm) {

	// first row
	data[0][0] = rm.cos_tet * rm.cos_psi;
	data[0][1] = rm.cos_tet * rm.sin_psi;
	data[0][2] = -1.0*rm.sin_tet;

	//second row
	data[1][0] = -1.0 * rm.sin_psi;
	data[1][1] = +1.0 * rm.cos_psi;
	data[1][2] = 0.0;

	// third row
	data[2][0] = rm.sin_tet * rm.cos_psi;
	data[2][1] = rm.sin_tet * rm.sin_psi;
	data[2][2] = rm.cos_tet;

};
// inverse rotation matrix 
// (transpose of direct since rot matrices are orthogonal)
void t_SqMat3::set_inv(const t_MatRotN& rm) {

	// first row
	data[0][0] = rm.cos_tet * rm.cos_psi;
	data[0][1] = -1.0 * rm.sin_psi;
	data[0][2] = rm.sin_tet * rm.cos_psi;

	//second row
	data[1][0] = rm.cos_tet * rm.sin_psi;
	data[1][1] = +1.0 * rm.cos_psi;
	data[1][2] = rm.sin_tet * rm.sin_psi;

	// third row
	data[2][0] = -1.0 * rm.sin_tet;
	data[2][1] = 0.0;
	data[2][2] = rm.cos_tet;

};

//  Description:                                                              //
//     This routine uses Crout's method to decompose a row interchanged       //
//     version of the n x n matrix A into a lower triangular matrix L and a   //
//     unit upper triangular matrix U such that A = LU.                       //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.                                                        //
//     Note!  In Crout's method the diagonal elements of U are 1 and are      //
//            not stored.                                                     //
//     Note!  The determinant of A is the product of the diagonal elements    //
//            of L.  (det A = det L * det U = det L).                         //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is                         //
//                      Crout_LU_with_Pivoting_Solve.                         //
//     (see below).                                                           //
//                                                                            //
//     The Crout method with partial pivoting is: Determine the pivot row and //
//     interchange the current row with the pivot row, then assuming that     //
//     row k is the current row, k = 0, ..., n - 1 evaluate in order the      //
//     the following pair of expressions                                      //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                                 for i = k, ... , n-1,                      //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                                                  / L[k][k] //
//                                      for j = k+1, ... , n-1.               //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A       Pointer to the first element of the matrix A[n][n].    //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     int     n       The number of rows or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //

int Crout_LU_Decomposition_with_Pivoting(double* A, int pivot[], int n)
{

	int row, i, j, k, p;
	double* p_k, * p_row, * p_col;
	double max;

	p_row = nullptr;
	p_col = nullptr;

	//         For each row and column, k = 0, ..., n-1,

	for (k = 0, p_k = A; k < n; p_k += n, k++) {

		//            find the pivot row

		pivot[k] = k;
		max = fabs(*(p_k + k));
		for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n) {
			if (max < fabs(*(p_row + k))) {
				max = fabs(*(p_row + k));
				pivot[k] = j;
				p_col = p_row;
			}
		}

		//     and if the pivot row differs from the current row, then
		//     interchange the two rows.

		if (pivot[k] != k)
			for (j = 0; j < n; j++) {
				max = *(p_k + j);
				*(p_k + j) = *(p_col + j);
				*(p_col + j) = max;
			}

		//                and if the matrix is singular, return error

		if (*(p_k + k) == 0.0) return -1;

		//      otherwise find the upper triangular matrix elements for row k. 

		for (j = k + 1; j < n; j++) {
			*(p_k + j) /= *(p_k + k);
		}

		//            update remaining matrix

		for (i = k + 1, p_row = p_k + n; i < n; p_row += n, i++)
			for (j = k + 1; j < n; j++)
				*(p_row + j) -= *(p_row + k) * *(p_k + j);

	}
	return 0;
}

int Crout_LU_with_Pivoting_Solve(double* LU, double B[], int pivot[],
	double x[], int n)
{
	int i, k;
	double* p_k;
	double dum;

	//         Solve the linear equation Lx = B for x, where L is a lower
	//         triangular matrix.                                      

	for (k = 0, p_k = LU; k < n; p_k += n, k++) {
		if (pivot[k] != k) { dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum; }
		x[k] = B[k];
		for (i = 0; i < k; i++) x[k] -= x[i] * *(p_k + i);
		x[k] /= *(p_k + k);
	}

	//         Solve the linear equation Ux = y, where y is the solution
	//         obtained above of Lx = B and U is an upper triangular matrix.
	//         The diagonal part of the upper triangular part of the matrix is
	//         assumed to be 1.0.

	for (k = n - 1, p_k = LU + n * (n - 1); k >= 0; k--, p_k -= n) {
		if (pivot[k] != k) { dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum; }
		for (i = k + 1; i < n; i++) x[k] -= x[i] * *(p_k + i);
		if (*(p_k + k) == 0.0) return -1;
	}

	return 0;
}

void test_LU() {
	t_SqMat3 mat;
	// simple test
	//for (int i = 0; i < 3; i++) mat[i][i] = i + 1;
	// more complex test with zeros on diagonal 
	mat[0][1] = 1.0;
	mat[1][2] = 1.0;
	mat[2][0] = 1.0;
	// insert your crazy test here )
	mat[0][0] = 1;
	mat[0][1] = 2;
	mat[0][2] = 3;

	mat[1][0] = 3;
	mat[1][1] = 4;
	mat[1][2] = 5;

	mat[2][0] = 7;
	mat[2][1] = 10;
	mat[2][2] = 0.02;

	hsLogMessage("MAt=%s", mat.to_str().c_str());
	t_SqMat3 mat_inv;
	mat_inv = mat.CalcInv();
	t_SqMat3 res;
	hsLogMessage("MAt_inv=%s", mat_inv.to_str().c_str());
	res = mat*mat_inv;
	hsLogMessage("Mat*MAt_inv=%s", res.to_str().c_str());

};

