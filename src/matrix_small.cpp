#include "matrix_small.h"

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