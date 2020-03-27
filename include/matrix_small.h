#pragma once

#include <array>

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
};

template<int N> t_Vec<N> operator*(double val, const t_Vec<N>& vec) {
	return vec*val;
};

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
