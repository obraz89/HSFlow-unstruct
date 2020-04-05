#pragma once

#include <map>

#include "cgnslib.h"

#include "common_data.h"
#include "common_procs.h"

using t_BufCGSize = t_BufInds<cgsize_t>;

// cgns section
struct t_CGSection {

	CG_ElementType_t itype;

	std::string name;

	// container to store indices of vertices
	t_BufCGSize _buf;
	// range of 1-based cgns indices for the section
	cgsize_t id_start, id_end;

	t_BufCGSize& get_buf() { return _buf; }
	const t_BufCGSize& get_buf() const{ return _buf; }

	cgsize_t* get_buf_data() { return _buf.data(); }
	const cgsize_t* get_buf_data() const{ return _buf.data(); }

	void alloc(cgsize_t n_elems, int n_verts_in_elem) {

		_buf.allocate(n_elems, n_verts_in_elem);
	}

	t_CGSection() :
		_buf(0, 0), itype(CG_ElementTypeNull), id_start(0), id_end(0), name("") {}
	t_CGSection(const char* sect_name) :
		_buf(0, 0), itype(CG_ElementTypeNull), id_start(0), id_end(0), name(sect_name) {}

	~t_CGSection() {  }


};

enum struct t_CGSectionKind{
	Cell=0,
	BC,
	Abutted,
	Undefined
};

class t_CGNSZone
{
	// vector of arrays of cells {[Tetra], [HEXA-8], ...}
	std::vector<t_CGSection*> pSectsCell;

	// bc patches
	std::vector<t_CGSection*> pSectsBC;
	// abutting patches
	std::vector<t_CGSection*> pSectsAbut;

	std::vector<t_CGSection*> pSectsAll;
public:
	const std::vector<t_CGSection*>& getSectsCell() const { return pSectsCell; }
	const std::vector<t_CGSection*>& getSectsBC() const { return pSectsBC; }
	const std::vector<t_CGSection*>& getSectsAbut() const { return pSectsAbut; }

	t_CGSection& getSectionCell(int i) { return *pSectsCell[i]; }
	const t_CGSection& getSectionCell(int i) const{ return *pSectsCell[i]; }

	t_CGSection& getSectionBC(int i) { return *pSectsBC[i]; }
	const t_CGSection& getSectionBC(int i) const { return *pSectsBC[i]; }

	t_CGSection& getSectionAbut(int i) { return *pSectsAbut[i]; }
	const t_CGSection& getSectionAbut(int i) const { return *pSectsAbut[i]; }


	cgsize_t countCells() {
		cgsize_t N=0;
		for (int i = 0; i < pSectsCell.size(); i++) 
			N += pSectsCell[i]->get_buf().nRows;
		return N;
	}

	void addSection(t_CGSection* psect, t_CGSectionKind kind) {

		switch (kind)
		{
		case t_CGSectionKind::Cell:
			pSectsCell.push_back(psect);
			break;
		case t_CGSectionKind::BC:
			pSectsBC.push_back(psect);
			break;
		case t_CGSectionKind::Abutted:
			pSectsAbut.push_back(psect);
			break;
		default:
			hsLogWarning("t_CGNSZone::addSection: adding unknown section...");
			break;
		}

		pSectsAll.push_back(psect);

	}

	t_CGNSZone() : pSectsCell(), pSectsBC(), pSectsAbut() { ; }
	~t_CGNSZone() { 
		for (int i = 0; i < pSectsCell.size(); i++) delete pSectsCell[i];
		for (int i = 0; i < pSectsAbut.size(); i++) delete pSectsAbut[i];
		for (int i = 0; i < pSectsBC.size(); i++) delete pSectsBC[i];
	}
};

struct t_CGNSContext
{
	int iFile, iBase;
	std::map<std::string, int>  map_ZneName2Idx;

	t_CGNSZone* cgZones;

	t_CGNSContext() : cgZones(NULL) { ; }
	~t_CGNSContext() { delete[] cgZones; }
};

int read_cgns_mesh();