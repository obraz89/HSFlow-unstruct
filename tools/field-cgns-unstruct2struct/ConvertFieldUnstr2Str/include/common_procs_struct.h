///////////////////////////////////////////////////////////////////////////////
// Name:        common_procs.h
// Purpose:     Generic procedures shared through main exe and plugins
// Author:      Andrey V. Novikov
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <ctime>
//-----------------------------------------------------------------------------

//
// HSFlow logging
//
#include "logging.h"
//-----------------------------------------------------------------------------

//
// String manipulation
//
std::string hs_string_format( const char* fmt, ... )
#if defined(__GNUC__)
	__attribute__(( format(printf, 1, 2) ))  // ask GCC to type-check against a format string
#endif
;

bool hs_string_ends_with(const std::string& fullString, const std::string& ending);
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//
// File system utilities
//

//-----------------------------------------------------------------------------

//
// Data arrays handling
//

/// Reimplementation std::copy_n() for older compilers not supporting it
template<typename T>
inline void copy_n(const T* from, size_t count, T* to)
{
	std::copy(from, from+count, to);
}


/**
 * Several arrays packed (stacked up) into one. Static version.
 * Each subarray consists of several POD elements.
 * Used e.g. for MPI trasfer
 *
 * @tparam T  type of the elements
 * @tparam N  number of subarrays
 * @tparam M  number of elements per subarray
 */
template<typename T, unsigned N, unsigned M>
class TpakArrays
{
	T _data[N * M];

public:
	TpakArrays() = default;
	TpakArrays(std::initializer_list<T> lst) {
		assert(lst.size() == N*M);
		std::copy(lst.begin(), lst.end(), _data);
	}

	T* data() {
		return _data;
	}
	/// Total elements count in all subarrays
	unsigned size() const {
		return N * M;
	}

	/// Subarray pointer
	T* operator[](const unsigned m) {
		return _data + m * M;
	}
	/// Subarray elements count
	unsigned subsize() const {
		return M;
	}
};


/**
 * Several arrays packed (stacked up) into one. Dynamic version.
 * Each subarray consists of several POD elements.
 *
 * @tparam T - type of the element, POD
 */
template<typename T>
class TpakArraysDyn
{
	T* _data = nullptr;
	unsigned N = 0;   // number of subarrays
	unsigned M = 0;   // number of elements per subarray

	TpakArraysDyn(TpakArraysDyn&) = delete;
	void operator=(TpakArraysDyn&) = delete;

public:
	TpakArraysDyn() = default;
	TpakArraysDyn(unsigned aN, unsigned aM) : N(aN), M(aM) {
		_data = new T[size()];
	}
	void reset(unsigned aN, unsigned aM) {
		delete[] _data;
		N = aN;  M = aM;
		_data = new T[size()];
	}

	T* data() {
		return _data;
	}
	/// Total elements count in all subarrays
	unsigned size() const {
		return N * M;
	}

	/// Subarray pointer
	T* operator[](const unsigned m) {
		return _data + m * M;
	}
	/// Subarray elements count
	unsigned subsize() const {
		return M;
	}
	~TpakArraysDyn(){
		delete[] _data;
	}
};
//-----------------------------------------------------------------------------

//
// Math functions
//

/// Squared value
template<typename T>
inline T sq(T x)
{
	return x*x;
}
//-----------------------------------------------------------------------------
