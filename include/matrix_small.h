#pragma once

#include <array>

#include <string>

#include <sstream>

#include "logging.h"

// classes for small matrix computations
// small means that size is known at compile time and is ... small

template<int N> class t_Vec {
protected:
	std::array<double,N> data;
public:
	t_Vec() { for (int i = 0; i < N; i++) data[i] = 0.0; }
	t_Vec(const double(&arr)[N]) { for (int i = 0; i < N; i++) data[i] = arr[i]; }

	double& operator[](int ind) { return data[ind]; }
	const double& operator[](int ind) const{ return data[ind]; }

	t_Vec& operator+=(const t_Vec& v) { for (int i = 0; i < N; i++) data[i] += v[i]; return *this; }
	t_Vec& operator*=(const double c) { for (int i = 0; i < N; i++) data[i] *= c; return *this; }

	t_Vec operator+(const t_Vec& v) const { t_Vec ret; for (int i = 0; i < N; i++) ret[i] = data[i] + v[i]; return ret; }
	t_Vec operator-(const t_Vec& v) const { t_Vec ret; for (int i = 0; i < N; i++) ret[i] = data[i] - v[i]; return ret; }
	t_Vec operator*(const double c) const { t_Vec ret; for (int i = 0; i < N; i++) ret[i] = c*data[i]; return ret; }

	double sq() const { double ret=0.0; for (int i = 0; i < N; i++) ret += data[i] * data[i]; return ret; }
	double norm() const { return sqrt(sq()); }

	void reset() { for (int i = 0; i < N; i++) data[i] = 0.0; }

	void cpy_to_c_arr(double(&c_arr)[N]) const {
		for (int i = 0; i < N; i++) c_arr[i] = data[i];
	}

	/**
	*  Make vector to be of unit length
	*  @return previous vector length
	*/
	void normalize() {double d = norm();for (int i = 0; i < N; i++) data[i] /= d;}
	void flip() { for (int i = 0; i < N; i++) data[i] *= -1; }

	/// Scalar product
	double dot(const t_Vec& v) const {
		double ret=0.0; 
		for (int i = 0; i < N; i++) 
			ret += data[i] * v[i]; 
		return ret;
	}

	std::string to_str() const{
		std::ostringstream ostr;
		ostr << "{";
			for (int i = 0; i < N; i++) ostr << data[i] << ";";
		ostr << "}\n";
		return ostr.str();
	}
};

template<int N> t_Vec<N> operator*(double val, const t_Vec<N>& vec) {
	return vec*val;
};

// this version of sign return 0 for 0 value
template <typename T> int sign(T val) {
	return (T(0) < val) - (val < T(0));
}

class t_Vec3 :public t_Vec<3> {
public:
	t_Vec3() :t_Vec<3>() {}
	t_Vec3(double a1, double a2, double a3) : t_Vec<3>({a1,a2,a3}) { ; }
	t_Vec3(const t_Vec<3>& v) :t_Vec<3>(v) {}

	void set(double a1, double a2, double a3) { data[0] = a1;  data[1] = a2;  data[2] = a3; }

	t_Vec3 cross(const t_Vec3& v) const	{
		return t_Vec3(
			data[1]*v[2] - data[2]*v[1], 
			data[2]*v[0] - data[0]*v[2], 
			data[0]*v[1] - data[1]*v[0]);	
	}

};

t_Vec3 operator*(double val, const t_Vec3& vec);

// rotation matrix in 3D space (can be used to rotate fluxes too)
// produced by a normal vector N
// the storage is minimized,
// trigonometrical calculations also minimize
// Convention:
// local Reference frame ked (ksi, eta, dzeta), global - XYZ
// coordinates of N in global RF: (Nx, Ny, Nz)
// in local RF: (Nk, Ne, Nd)
// unity vectors of axis of ked in global RF: (ek, ee, ed)
// The vector N generates x-axis of local RF : ek = N
// ee and ed is an arbitrary pair of vectors normal to ek.
// We choose them with Psi and Theta angles (Yaw and Pitch)
// We define direct rotation matrix R as follows (_T is transposition)
// (k,e,d)_T = R*(X,Y,Z)_T
struct t_MatRotN {

	double sin_psi, cos_psi;
	double sin_tet, cos_tet;

	void calc_rot_angles_by_N(const t_Vec3& N);

};



int Crout_LU_Decomposition_with_Pivoting(double* A, int pivot[], int n);

int Crout_LU_with_Pivoting_Solve(double* LU, double B[], int pivot[],
	double x[], int n);

// matrix
// packaged by rows, A[1][2] - second row, third column
template<int NRows, int NCols> class t_Mat {
protected:
	std::array<std::array<double, NCols>, NRows> data;

	void _mul_by_vec(const t_Vec<NCols>& vec, t_Vec<NRows>& dest) const {
		for (int i = 0; i < NRows; i++) {
			dest[i] = 0.0;
			for (int j = 0; j < NCols; j++) dest[i] += data[i][j] * vec[j];
		}
	};

	template<int NColsRgt>
	void _mul_by_mat(const t_Mat<NCols, NColsRgt>& rMat, t_Mat<NRows, NColsRgt>& dest) const{
		dest.reset();
		for (int i = 0; i < NRows; i++)
			for (int j = 0; j < NColsRgt; j++)
				for (int k = 0; k < NCols; k++)
					dest[i][j] += data[i][k] * rMat[k][j];
	}

public:
	t_Mat() { 
		for (int i = 0; i < NRows; i++) 
			for (int j = 0; j < NCols; j++) 
				data[i][j] = 0.0; 
	};
	void add(const t_Mat& dM) {
		for (int i = 0; i < NRows; i++)
			for (int j = 0; j < NCols; j++)
				data[i][j] += dM[i][j];
	}
	std::array<double, NCols>& operator[](int iRow) {
		return data[iRow];
	};

	const std::array<double, NCols>& operator[](int iRow) const {
		return data[iRow];
	};

	t_Vec<NRows>  operator*(const t_Vec<NCols>& vec) const {
		t_Vec<NRows> ret;
		_mul_by_vec(vec, ret);
		return ret;
	}

	void setRow(int iRow, const t_Vec<NCols>& vec) {

		if ((iRow<0) || (iRow>NRows - 1)) 
			hsLogMessage("Mat:setRow:Error: index out of range: row=%d", iRow);

		for (int j = 0; j < NCols; j++) data[iRow][j] = vec[j];

	}

	void setCol(int iCol, const t_Vec<NRows>& vec) {

		if ((iCol<0) || (iCol>NCols - 1))
			hsLogMessage("Mat:setRow:Error: index out of range: row=%d", iCol);

		for (int i = 0; i < NRows; i++) data[i][iCol] = vec[i];

	}
	// primitive l2 norm
	// square root of sum of squares of all elements
	double calcNormPrimitive() {

		double norm = 0.0;

		for (int i = 0; i < NRows; i++)
			for (int j = 0; j < NCols; j++)
				norm += data[i][j]*data[i][j];
		return sqrt(norm);

	}

	int getFlatInd(int i, int j) const{ return i * NCols + j; }

	void flatten(double(&p)[NRows * NCols]) const {
		for (int i = 0; i < NRows; i++)
			for (int j = 0; j < NCols; j++)
				p[getFlatInd(i, j)] = data[i][j];
	}

	std::string to_str() {
		std::ostringstream ostr;
		ostr << "Mat:";
		for (int i = 0; i < NRows; i++) {
			ostr << "[";
			for (int j = 0; j < NCols; j++)
				ostr << data[i][j] << ";";
			ostr << "]\n";
		}
		return ostr.str();
	}

	void reset() {
		for (int i = 0; i < NRows; i++)
			for (int j = 0; j < NCols; j++)
				data[i][j] = 0.0;
	}
};

template<int N> class t_SqMat : public t_Mat<N,N> {
public:
	t_SqMat():t_Mat<N,N>(){};

	void setToUnity() {
		t_Mat<N, N>::reset();
		for (int i = 0; i < N; i++) t_Mat<N, N>::data[i][i] = 1.0;
	}

	void LU(t_SqMat& L, t_SqMat& U) {

		double A[N * N];

		flatten(A);

		int pivot[N];

		Crout_LU_Decomposition_with_Pivoting(A, pivot, N);

		L.reset();

		for (int i = 0; i < N; i++)
			for (int j = 0; j <= i; j++)
				L[i][j] = A[t_Mat<N, N>::getFlatInd(i,j)];

		U.setToUnity();

		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
				U[i][j] = A[t_Mat<N, N>::getFlatInd(i, j)];

	}

	t_SqMat CalcInv() const{

		double A[N * N];

		int pivot[N];

		t_Mat<N,N>::flatten(A);

		int ok = Crout_LU_Decomposition_with_Pivoting(A, pivot, N);

		if (ok < 0)
			hsLogMessage("Error: t_SqMat::CalcInv: matrix is singular");

		t_SqMat inv;

		double B[N];

		double x[N];

		for (int j = 0; j < N; j++) {

			for (int i = 0; i < N; i++) B[i] = 0.0;
			B[j] = 1.0;

			Crout_LU_with_Pivoting_Solve(A, B, pivot, x, N);

			for (int i = 0; i < N; i++) inv[i][j] = x[i];

		}

		return inv;

	}

	t_SqMat<N> operator*(const t_SqMat<N>& rMat) const {
		t_SqMat<N> ret;
		t_Mat<N, N>::_mul_by_mat<N>((const t_Mat<N, N>&)rMat, (t_Mat<N, N>&)ret);
		return ret;
	}

	t_Vec<N> operator*(const t_Vec<N>& rVec) const {
		t_Vec<N> ret;
		t_Mat<N,N>::_mul_by_vec((const t_Vec<N>&)rVec, (t_Vec<N>&)ret);
		return ret;
	}

};

class t_SqMat3 :public t_SqMat<3> {
public:
	t_SqMat3() :t_SqMat<3>() {}
	void set(const t_MatRotN& rot_mat);
	void set_inv(const t_MatRotN& rot_mat);
	void operator=(const t_SqMat<3>& a_mat) { (t_SqMat<3>&)* this = a_mat; }
	t_Vec3 operator*(const t_Vec3& v) const{
		t_Vec3 ret;
		_mul_by_vec(v, ret);
		return ret;
	}

	t_SqMat3 operator*(const t_SqMat3& rmat) const {
		t_SqMat3 ret;
		_mul_by_mat<3>((const t_Mat<3,3>&)rmat, (t_Mat<3,3>&)ret);
		return ret;

	}
};

void test_LU();
