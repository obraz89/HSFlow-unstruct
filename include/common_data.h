#pragma once

#include <map>

#include <array>
#include <vector>

#include "logging.h"

#include "matrix_small.h"

typedef __int64 lint;

enum struct t_FaceBCKind {
	Fluid = 0,
	Inflow,
	Outflow,
	Sym,
	Wall
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

//
// Problem solving state
//
struct TState
{
	int mpiRank, mpiNProcs;  // MPI rank, number of procs
							 // Relation of zone index and MPI rank working with it
	int* map_zone2rank;  // map_zone2rank[izne] == mpiRank, where izne -- 0-based zone index

						 //---

	int nTmStep;     // current time step number
	int nwtIter;     // current Newton's iteration number

					 /// L_inf norm of residual at 0-th Newton's iteration at current step
	double initialResidual;
};

extern TState G_State;
