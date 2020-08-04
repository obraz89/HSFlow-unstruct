///////////////////////////////////////////////////////////////////////////////
// Name:        evector.h
// Purpose:     Euclidean (geometric) vectors, 2D & 3D
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>
//-----------------------------------------------------------------------------

/**
 * Euclidean vector dummy template
 */
template<bool is3D>
struct Tvec;
//-----------------------------------------------------------------------------


/**
 * Euclidean 2D vector
 */
template<>
struct Tvec<false>
{
	double a, b;

	Tvec() = default; // zero initialization
	Tvec(double aa, double ab) : a(aa), b(ab) { ; }

	void set(double aa, double ab){ a = aa;  b = ab; }

	void operator+=(const Tvec& v){ a += v.a;  b += v.b;  }
	Tvec operator+(const Tvec& v) const { return Tvec{a + v.a, b + v.b}; }
	Tvec operator-(const Tvec& v) const { return Tvec{a - v.a, b - v.b}; }
	void operator*=(double k){ a *= k;  b *= k;  }
	Tvec operator*(double k) const { return Tvec{a * k, b * k}; }

	double sq() const {  return a*a + b*b;  }  //< Square
	double length() const {  return sqrt( sq() );  }

	/**
	 *  Make vector to be of unit length
	 *  @return previous vector length
	 */
	double normalize() {
		const double d = length();
		a /= d;  b /= d;
		return d;
	}

	void flip(){  a = -a;  b = -b;  }

	Tvec orthogonal() const {
		return Tvec(-b, a);
	}

	/// Scalar product
	double dot(const Tvec& v) const {
		return a*v.a + b*v.b;
	}
};
//-----------------------------------------------------------------------------


/**
 * Euclidean 3D vector
 */
template<>
struct Tvec<true>
{
	double a, b, c;

	Tvec() = default; // zero initialization
	Tvec(double a1, double a2, double a3) : a(a1), b(a2), c(a3) { ; }

	void set(double a1, double a2, double a3){ a = a1;  b = a2;  c = a3; }

	void operator+=(const Tvec& v){ a += v.a;  b += v.b;  c += v.c; }
	Tvec operator+(const Tvec& v) const { return Tvec{a + v.a, b + v.b, c + v.c}; }
	Tvec operator-(const Tvec& v) const { return Tvec{a - v.a, b - v.b, c - v.c}; }
	void operator*=(double k){ a *= k;  b *= k;  c *= k; }
	Tvec operator*(double k) const { return Tvec{a * k, b * k, c * k}; }

	double sq() const {  return a*a + b*b + c*c;  }
	double length() const {  return sqrt( sq() );  }

	/**
	 *  Make vector to be of unit length
	 *  @return previous vector length
	 */
	double normalize() {
		const double d = length();
		a /= d; b /= d; c /= d;
		return d;
	}
	void flip() {  a = -a;  b = -b;  c = -c;  }

	/// Vector product
	Tvec cross(const Tvec& v) const
	{
		return Tvec(b*v.c - c*v.b, c*v.a - a*v.c, a*v.b - b*v.a);
	}

	/// Scalar product
	double dot(const Tvec& v) const
	{
		return a*v.a + b*v.b + c*v.c;
	}
};
//-----------------------------------------------------------------------------

template<bool is3D>
Tvec<is3D> operator*(double k, Tvec<is3D> v) { v *= k; return v; }

/// Euclidean vector aliases
using Tvec3D = Tvec<true>;
using Tvec2D = Tvec<false>;
//-----------------------------------------------------------------------------
