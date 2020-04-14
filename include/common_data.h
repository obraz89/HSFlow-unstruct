#pragma once

#include <map>

#include <array>
#include <vector>

#include "logging.h"

#include "matrix_small.h"

typedef __int64 lint;

// identifier of the particular bc reduced to int
// usually bc ids count as enum 0,1,2,3,4
// using magic value -7 for id of fluid cell not to intersect with bc sets
class t_FaceBCID {
	int val;
public:
	static const int Fluid = -7;
	void set(int a_val) { val = a_val; }
	int get() const { return val; }
	bool operator==(int vv) const{ return val == vv; }

};

enum struct t_CellKind {
	Tetra = 0,
	Brick,
	Pyra,
	None
};

template<typename T>
class t_BufInds {
	T* buf;
public:
	T nRows, nCols;
	T* data() { return buf; }
	const T* data() const { return buf; }
	t_BufInds() = delete;
	t_BufInds(t_BufInds&) = delete;
	t_BufInds(T a_NR, T a_NC) { nRows = a_NR; nCols = a_NC; buf = new T[nRows*nCols]; }
	void allocate(T a_NR, T a_NC) { delete[] buf;  nRows = a_NR; nCols = a_NC; buf = new T[nRows*nCols]; };
	T& get_val(T i, T j) { return *(buf + i*nCols + j); };
	const T& get_val(T i, T j) const { return *(buf + i*nCols + j); }
	~t_BufInds() { delete[] buf; }

};

// plain set of indices, usually decomposition of a face
template<typename T, int Nmax>
class t_TSet {
	std::array<T, Nmax> buf;
	int NElems;
public:
	t_TSet() :NElems(0) { for (int i = 0; i < Nmax; i++) buf[i] = 0; }
	int size() const { return NElems; }
	void setSize(int a_size) { 
#ifdef _DEBUG
		if (a_size > Nmax) hsLogMessage("t_TSet:setSize:Error: size is too big");
#endif // _DEBUG
		NElems = a_size;
	}
	T& operator[](int i) { 
#ifdef _DEBUG
		if ((i > Nmax-1) || (i<0)) hsLogMessage("t_TSet: wrong index");
#endif // _DEBUG
		return buf[i]; }
	const T& operator[](int i) const { 
#ifdef _DEBUG
		if ((i > Nmax - 1) || (i<0)) hsLogMessage("t_TSet: wrong index");
#endif // _DEBUG
		return buf[i]; }

	// weak comparison, order not important {1,2,3,4}=={1,3,4,2} : true
	static bool cmp_weak(const t_TSet<T, Nmax>& lv, const t_TSet<T, Nmax>& rv){
		if (lv.size() != rv.size()) return false;
		std::array<T, Nmax> l_sorted = lv.buf;
		std::sort(l_sorted.begin(), l_sorted.end());
		std::array<T, Nmax> r_sorted = rv.buf;
		std::sort(r_sorted.begin(), r_sorted.end());
		bool cmp_ok = true;
		for (int i = 0; i < Nmax; i++) cmp_ok = cmp_ok && (l_sorted[i] == r_sorted[i]);
		return cmp_ok;
	};

	static bool cmp_strict(const t_TSet<T, Nmax>& lv, const t_TSet<T, Nmax>& rv){
		if (lv.size() != rv.size()) return false;
		bool cmp_ok = true;
		for (int i = 0; i < NMax; i++) cmp_ok = cmp_ok && (l_sorted[i] == r_sorted[i]);
		return cmp_ok;
	}
};

// container for small packs of index sets
template<typename T, int nRows, int nCols>
class t_BufIndsStat {
private:
	T buf[nRows][nCols];
public:
	T* data() { return buf; }
	T& get_val(int i, int j) { return buf[i][j]; };
	const T& get_val(int i, int j) const { return buf[i][j]; }
};

/**
 * Double array wrapper
 */
class t_ArrDbl
{
	double* _D;
	size_t _size;

	t_ArrDbl(t_ArrDbl&) = delete;
	void operator=(t_ArrDbl&) = delete;

public:
	t_ArrDbl() : _D(nullptr), _size(0) { ; }

	void alloc(size_t size) {
		_D = new double[size];
		_size = size;
	}

	size_t size() const { return _size; }

	void set(size_t n, const double& val) {
		_D[n] = val;
	}

	double* data() {
		return _D;
	}

	//MPI_Datatype typeMPI() const {
	//	return MPI_DOUBLE;
	//}

	~t_ArrDbl() {
		delete[] _D;
	}
};

//
// Problem solving state
//
struct TState
{
	int mpiRank, mpiNProcs;  // MPI rank, number of procs
							 // Relation of zone index and MPI rank working with it
	int* map_zone2rank;  // map_zone2rank[izne] == mpiRank, where izne -- 0-based zone index

						 //---

	//int nTmStep;     // current time step number
	//int nwtIter;     // current Newton's iteration number

					 /// L_inf norm of residual at 0-th Newton's iteration at current step
	//double initialResidual;

	// TODO: these are scheme parameters, separate them later

	double time;

	double ResidTot;

	// debug vars
	double ResidNormVeloWall;
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
	~TpakArraysDyn() {
		delete[] _data;
	}
};

extern TState G_State;
