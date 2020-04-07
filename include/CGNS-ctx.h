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

	t_CGSection() = delete;
	t_CGSection(const char* sect_name, cgsize_t istart, cgsize_t iend) :
		_buf(0, 0), itype(CG_ElementTypeNull), id_start(istart), id_end(iend), name(sect_name) {}

	~t_CGSection() {  }


};

enum struct t_CGSectionKind{
	Cell=0,
	BC,
	Abutted,
	Undefined
};


// Element connectivity set
// 1-to-1 connectivity
// IMPORTANT: this is NOT vertex connectivity (used in structured code)
// it is more general element-to-element connectivity
// _buf stores pairs of connections: {id_my, id_dnr}
class t_CGConnSet {
	cgsize_t idZoneMy, idZoneDnr;
	t_BufCGSize _buf;

public:
	t_BufCGSize& get_buf() { return _buf; }
	const t_BufCGSize& get_buf() const { return _buf; }

	cgsize_t* get_buf_data() { return _buf.data(); }
	const cgsize_t* get_buf_data() const { return _buf.data(); }

	cgsize_t getZoneIDMy() const{ return idZoneMy; }
	cgsize_t getZoneIDDnr() const{ return idZoneDnr; }

	void getConnIds(cgsize_t i, cgsize_t& id_my, cgsize_t& id_dnr) const{
		id_my = _buf.get_val(0, i);
		id_dnr = _buf.get_val(1, i);
	};

	t_CGConnSet() = delete;
	t_CGConnSet(cgsize_t id_Z_My, cgsize_t id_Z_Dnr, cgsize_t N):
		idZoneMy(id_Z_My), idZoneDnr(id_Z_Dnr), _buf(0,0) {
		_buf.allocate(2, N);
	}

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

	std::vector<t_CGConnSet*> pConns;

	std::string name;

	int NVerts, NCells;

	double* x_coords, * y_coords, * z_coords;

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

	const std::vector<t_CGConnSet*>& getConns() const { return pConns; }

	void setNVertsNCells(cgsize_t nv, cgsize_t nc) { NVerts = nv; NCells = nc; }
	cgsize_t getNVerts() const { return NVerts; }
	cgsize_t getNCells() const { return NCells; }

	double*& getXCoords() { return x_coords; }
	double*& getYCoords() { return y_coords; }
	double*& getZCoords() { return z_coords; }

	void setName(const char* a_name) { name = std::string(a_name); };
	const std::string& getName() const{ return name; };


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

	std::vector<cgsize_t> getVertsOfElem(cgsize_t id) const{

		std::vector<cgsize_t> vert_ids;

		for (int i = 0; i < pSectsAll.size(); i++) {

			const t_CGSection& sect = *pSectsAll[i];

			if (sect.id_start <= id && id <= sect.id_end) {

				int nverts = sect.get_buf().nCols;

				vert_ids.resize(sect.get_buf().nCols);

				cgsize_t irow = id - sect.id_start;

				for (int j = 0; j < nverts; j++)
					vert_ids[j] = sect.get_buf().get_val(irow, j);

				return vert_ids;

			}
		}

		hsLogError("t_CGNSZone:getVertsOfElem: can't find the element with id=%ld", id);
		return vert_ids;

	};

	void addConn(t_CGConnSet* pCon) { pConns.push_back(pCon); }

	t_CGNSZone() : pSectsCell(), pSectsBC(), pSectsAbut(), name("") { ; }
	~t_CGNSZone() { 
		for (int i = 0; i < pSectsCell.size(); i++) delete pSectsCell[i];
		for (int i = 0; i < pSectsAbut.size(); i++) delete pSectsAbut[i];
		for (int i = 0; i < pSectsBC.size(); i++) delete pSectsBC[i];
		for (int i = 0; i < pConns.size(); i++) delete pConns[i];

		delete[] x_coords, y_coords, z_coords;
	}
};

struct t_CGNSContext
{
	int iFile, iBase;
	std::map<std::string, int>  map_ZneName2Idx;

	int nZones;

	t_CGNSZone* cgZones;

	t_CGNSContext() : cgZones(NULL) { ; }
	~t_CGNSContext() { delete[] cgZones; }

	cgsize_t getCGZoneIDByName(const char* name) {
		return map_ZneName2Idx[name] + 1;
	};

	cgsize_t getNumOfGhostsForZone(int cgZoneID) const;

	bool _parseConnectivity();

	bool checkBCs();

	bool loadGridCoords();

	bool readMesh(std::string fn);

	bool writeMesh();

	bool readField();

	bool writeField();
	

};

extern t_CGNSContext G_CGNSCtx;