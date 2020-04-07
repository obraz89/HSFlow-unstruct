#pragma once

#include <array>

#include <string>

#include <sstream>

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

	std::string to_str() {
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

// Square matrix
// packaged by rows, A[1][2] - second row, third column
template<int N> class t_SqMat {
protected:
	std::array<std::array<double,N>, N> data;

	void _mul_by_vec(const t_Vec<N>& vec, t_Vec<N>& dest) const{
		for (int i = 0; i < N; i++) {
			dest[i] = 0.0;
			for (int j = 0; j < N; j++) dest[i] += data[i][j] * vec[j];
		}
	};
public:
	t_SqMat(){ for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) data[i][j] = 0.0; };
	std::array<double, N>& operator[](int row) { return data[row]; }
	const std::array<double, N>& operator[](int row) const{ return data[row]; }


	t_Vec<N>  operator*(const t_Vec<N>& vec) const{
		t_Vec<N> ret;
		_mul_by_vec(vec, ret);
		return ret;
	}

	// debugging
	std::string to_str() {
		std::ostringstream ostr;
		ostr << "R:";
		for (int i = 0; i < N; i++) {
			ostr << "[";
			for (int j = 0; j < N; j++)
				ostr << data[i][j] << ";";
			ostr << "]\n";
		}
		return ostr.str();
	}

};

class t_SqMat3 :public t_SqMat<3> {
public:
	t_SqMat3() :t_SqMat<3>() {}
	void set(const t_MatRotN& rot_mat);
	void set_inv(const t_MatRotN& rot_mat);
	t_Vec3 operator*(const t_Vec3& v) const{
		t_Vec3 ret;
		_mul_by_vec(v, ret);
		return ret;
	}
};
