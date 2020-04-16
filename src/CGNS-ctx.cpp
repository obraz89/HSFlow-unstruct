#include "CGNS-ctx.h"

#include "logging.h"

#include "bc_common.h"

t_CGNSContext G_CGNSCtx;

bool t_CGNSContext::readMesh(std::string gridFN) {

		cgsize_t isize[3];

		cgsize_t irmin, irmax, istart, iend;

		const int CG_MAX_NAME_LENGTH = 128 + 1/*terminating 0*/;

		char zonename[CG_MAX_NAME_LENGTH], sectionname[CG_MAX_NAME_LENGTH];

		CGNS_ENUMT(ElementType_t) itype;

		int res = CG_OK;

		char szName[CG_MAX_NAME_LENGTH];  // names in CGNS file

		t_CGNSContext ctx;

		if (cg_open(gridFN.c_str(), CG_MODE_READ, &this->iFile) != CG_OK)
		{
			hsLogMessage("Can't open grid file '%s' for reading (%s)",
				gridFN, cg_get_error());
		}

		// assume only one base
		this->iBase = 1; 

		//
		// Space dimensions
		//
		int dimCell = 0, dimPhys = 0;
		cg_base_read(this->iFile, this->iBase, szName, &dimCell, &dimPhys);
		if (dimCell != 3)
		{
			hsLogMessage("The grid is not for 3D problems");
		}

		//
		// Number of zones
		//
		cg_nzones(this->iFile, this->iBase, &this->nZones);

		// Temporary zones data used on CGNS file parsing
		this->cgZones = new t_CGNSZone[this->nZones];


		for (int iZone = 0; iZone < this->nZones; iZone++) {

			const int& cgZneID = iZone + 1;	
			t_CGNSZone& cgZne = this->cgZones[iZone];

			hsLogMessage("Reading zone#%d", iZone);

			CG_ZoneType_t type;  cg_zone_type(this->iFile, this->iBase, cgZneID, &type);
			bool isOk = (type == CG_Unstructured) ? 1 : 0;
			if (!isOk) hsLogMessage("Only unstructured grids are supported");

			cg_zone_read(this->iFile, this->iBase, cgZneID, zonename, isize);

			// CGNS documentation: midlevel/structural.html#zone
			// isize = {NVertex, NCell3D, NBoundVertex}
			cgZne.setName(zonename); 

			cgZne.setNVertsNCells(isize[0], isize[1]);

			const cgsize_t& nVerts = cgZne.getNVerts();
			const cgsize_t& nCells = cgZne.getNCells();

			hsLogMessage("Number of Verts:%d", nVerts);

			hsLogMessage("Number of Cells:%d", nCells);

			//hsLogMessage("Number of BC Verts:%d", isize[2]);

			this->map_ZneName2Idx[cgZne.getName()] = iZone;

			// reading sections

			int nsections, index_sect, nbndry, iparent_flag;
			cgsize_t iparentdata;

			cg_nsections(this->iFile, this->iBase, cgZneID, &nsections);

			hsLogMessage("Number of sections:%d", nsections);

			for (index_sect = 1; index_sect <= nsections; index_sect++)
			{
				cgsize_t n_elems, n_verts_in_elem = 0;

				cg_section_read(this->iFile, this->iBase, cgZneID, index_sect, sectionname,
					&itype, &istart, &iend, &nbndry, &iparent_flag);
				printf("\nReading section info...\n");
				printf("   section name=%s\n", sectionname);
				printf("   section type=%s\n", ElementTypeName[itype]);
				printf("   istart,iend=%i, %i\n", (int)istart, (int)iend);

				// TODO: universal way to detect sections containing cells
				// for now detect by type of elements
				//if (itype != CG_HEXA_8 && itype != CG_TETRA_4) continue;

				if (itype == CG_HEXA_8) n_verts_in_elem = 8;
				if (itype == CG_TETRA_4) n_verts_in_elem = 4;
				if (itype == CG_PENTA_6) n_verts_in_elem = 6;

				n_elems = iend - istart + 1;

				// bc
				if (itype == CG_QUAD_4) n_verts_in_elem = 4;
				if (itype == CG_TRI_3) n_verts_in_elem = 3;

				// reading face patches
				if (itype == CG_QUAD_4 || itype == CG_TRI_3) {

					t_CGSection* pPatch = new t_CGSection(sectionname, istart, iend);

					if (G_pBCList->has(sectionname) == true) {
						t_FaceBCID faceBCID = G_pBCList->getBCID(sectionname);
						pPatch->BCId = faceBCID;
						cgZne.addSection(pPatch, t_CGSectionKind::BC);
					}
					else
						cgZne.addSection(pPatch, t_CGSectionKind::Abutted);

					pPatch->alloc(n_elems, n_verts_in_elem);

					hsLogMessage("   reading face patch : section %d, Type: %s\n",
						index_sect, ElementTypeName[itype]);

					cg_elements_read(this->iFile, this->iBase, cgZneID, index_sect, pPatch->get_buf_data(), \
						& iparentdata);

					pPatch->itype = itype;

					// debug output of section
					//for (int i = 0; i < n_elems; i++) {
					//	for (int j = 0; j < n_verts_in_elem; j++)
					//		std::cout << pPatch->get_buf().get_val(i, j) << ";";
					//	std::cout << std::endl;
					//}


					continue;
				}
				// reading elements
				if (itype == CG_HEXA_8 || 
					itype == CG_TETRA_4 || 
					itype == CG_PENTA_6) {

					t_CGSection* pNewCellSet = new t_CGSection(sectionname, istart, iend);

					cgZne.addSection(pNewCellSet, t_CGSectionKind::Cell);

					pNewCellSet->alloc(n_elems, n_verts_in_elem);

					hsLogMessage("   reading elements : section %d, Type: %s\n",
						index_sect, ElementTypeName[itype]);

					cg_elements_read(this->iFile, this->iBase, cgZneID, index_sect, pNewCellSet->get_buf_data(), \
						& iparentdata);

					pNewCellSet->itype = itype;

					// debug output of section
					//for (int i = 0; i < n_elems; i++) {
					//	for (int j = 0; j < n_verts_in_elem; j++)
					//		std::cout << pNewCellSet->get_buf().get_val(i, j) << ";";
					//	std::cout << std::endl;
					//}
					continue;
				}

				hsLogMessage("Error: read_cgns_mesh(): unsupported section type");

			}

		}

		loadGridCoords();

		// update ctx with connectivity info
		if (!_parseConnectivity())
			return false;

		// checks with some cgns bc-specific funcs
		if (!checkBCs())
			return false;

		return 0;


}

bool t_CGNSContext::_parseConnectivity() {

	for (int iZne = 0; iZne < nZones; ++iZne)
	{
		const int& cgZneID = iZne + 1;
		t_CGNSZone& cgZne = cgZones[iZne];

		int n1to1 = 0;

		if (G_State.mpiRank == 0)
			cg_n1to1(this->iFile, this->iBase, cgZneID, &n1to1);

		int NPatchesAbut = 0;

		if (G_State.mpiRank == 0)
			cg_nconns(this->iFile, this->iBase, cgZneID, &NPatchesAbut);

		//ier = cg_nconns(int fn, int B, int Z, int* nconns)
		//cg_conn_info();
		//cg_conn_read();
		
		//MPI_Bcast(&cgZne.nPatches, 1, MPI_INT, 0, PETSC_COMM_WORLD);

		if (NPatchesAbut + n1to1 < 1 && nZones > 1) {
			hsLogError("Missing inter-zone 1-to-1 connectivity in '%s'#%d", cgZne.getName(), cgZneID);
			//return false;
		}

		if (n1to1 > 0 )
			hsLogError(
				"parseConnectivity: 1-to-1 vertex connectivity detected, \
				 remake grid with generalized connectivity (via faces' connectivity)");

		if (NPatchesAbut != cgZne.getSectsAbut().size()) {
			hsLogError("parseConnectivity:Wrong number of abutted patches");
			//return false;
		}

		// Loop through connectivity patches
		for (int cgPatchId = 1; cgPatchId <= NPatchesAbut; ++cgPatchId)
		{
			t_CGSection& cgFacePatch = cgZne.getSectionAbut(cgPatchId - 1);

			char connectname[128];
			CG_GridLocation_t location;
			CG_GridConnectivityType_t connect_type;
			CG_PointSetType_t ptset_type;
			cgsize_t npnts;

			char donorname[128];
			CG_ZoneType_t donor_zonetype;
			CG_PointSetType_t donor_ptset_type;
			CG_DataType_t donor_datatype;
			cgsize_t ndata_donor[128];
			int ier = cg_conn_info(this->iFile, this->iBase, cgZneID, cgPatchId, connectname,
				&location, &connect_type,
				&ptset_type, &npnts, donorname,
				&donor_zonetype, &donor_ptset_type,
				&donor_datatype, ndata_donor);

			cgsize_t cgZneDnrID = this->getCGZoneIDByName(donorname);

			t_CGConnSet* pConnNew = new t_CGConnSet(cgZneID, cgZneDnrID, npnts);

			// a little bit dirty reading:
			// read points and donor points into the same buf array but with
			// corresponding shifts
			cgsize_t* pnts = pConnNew->get_buf_data();
			cgsize_t* pnts_dnr = pnts + npnts;

			ier = cg_conn_read(this->iFile, this->iBase, cgZneID, cgPatchId,
				pnts, donor_datatype, pnts_dnr);

			cgZne.addConn(pConnNew);

			// debug 
			//cgsize_t id_my, id_dnr;
			//for (int i = 0; i < npnts; i++) {
			//	const t_CGConnSet& conn = *cgZne.getConns().back();
			//	conn.getConnIds(i, id_my, id_dnr);
			//	hsLogMessage("Connectivity elem:%ld <-> %ld", id_my, id_dnr);
			//}

			
		}

	}
	return true;
}

bool t_CGNSContext::checkBCs() { 

	// bcs
	int nBCs;

	const int CG_MAX_NAME_LENGTH = 128 + 1/*terminating 0*/;

	for (int iZne = 0; iZne < nZones; ++iZne)
	{
		const int& cgZneID = iZne + 1;
		t_CGNSZone& cgZne = this->cgZones[iZne];

		cg_nbocos(this->iFile, this->iBase, cgZneID, &nBCs);

		hsLogMessage("BC check .................Number of BC sets:%d", nBCs);

		for (int iBC = 1; iBC <= nBCs; ++iBC) {

			char szPatchName[CG_MAX_NAME_LENGTH];

			CG_BCType_t iBCtype;

			CG_PointSetType_t pntSetType;
			cgsize_t nPnts = 0; // number of points defining the BC region

								// Normals to the BC patch - UNUSED
			int iNorm[3]; // normal as index vector (computational coords)
			cgsize_t normListSize = 0;  CG_DataType_t normDataType; // normals in phys coords

			int nDatasets = 0; // number of datasets with additional info for the BC

			cg_boco_info(this->iFile, this->iBase, cgZneID, iBC,
				szPatchName, &iBCtype,
				&pntSetType, &nPnts,
				iNorm, &normListSize, &normDataType,
				&nDatasets
			);
			if (pntSetType != CG_PointRange && nPnts != 2)
			{
				hsLogMessage(
					"Boundary condition patch '%s'(#%d) of zone '%s'(#%d) isn't defined as point range",
					szPatchName, iBC, cgZne.getName(), cgZneID);
				szPatchName[0] = 0x3;  // 'end of text' code -> error indicator
				return false;
			}
			cgsize_t idxRng[2];
			cg_boco_read(this->iFile, this->iBase, cgZneID, iBC, idxRng, nullptr);

			hsLogMessage("BC with name %s has idxrng_start=%d and idxrng_end=%d", szPatchName, idxRng[0], idxRng[1]);

			char szBC[CG_MAX_NAME_LENGTH] = "";

			cg_goto(this->iFile, this->iBase, "Zone_t", cgZneID, "ZoneBC", 0, "BC_t", iBC, NULL);
			if (cg_famname_read(szBC) == CG_OK)
			{
				// Read "Fam_Descr_Name" generated by Pointwise 16.03
				if (cg_goto(this->iFile, this->iBase, szBC, 0, NULL) == CG_OK)
				{
					char szDescrName[CG_MAX_NAME_LENGTH] = "";
					char* szDescrText = nullptr;
					if (cg_descriptor_read(1, szDescrName, &szDescrText) == CG_OK)
					{
						if (strcmp(szDescrName, "Fam_Descr_Name") == 0)
						{
							strcpy_s(szBC, szDescrText);
							cg_free(szDescrText);
						}
					}
				}
			}
			else
			{
				strcpy_s(szBC, szPatchName);
			}
		}
		hsLogMessage("BC check .................Ok");
		// ~bcs

	}

	return true; 
}


static bool parseVCs(t_CGNSContext& ctx) { return true; }

bool t_CGNSContext::loadGridCoords() {

	for (int iZne = 0; iZne < nZones; ++iZne) {

		t_CGNSZone& cgZne = cgZones[iZne];

		const cgsize_t& nVerts = cgZne.getNVerts();

		const int cgZneID = iZne + 1;
		const char* coordNames[] = { "CoordinateX", "CoordinateY", "CoordinateZ" };

		cgZne.getXCoords() = new double[nVerts];
		cgZne.getYCoords() = new double[nVerts];
		cgZne.getZCoords() = new double[nVerts];

		cgsize_t irmin = 1; cgsize_t irmax = nVerts;

		cg_coord_read(this->iFile, this->iBase, cgZneID, "CoordinateX",
			CG_RealDouble, &irmin, &irmax, cgZne.getXCoords());
		cg_coord_read(this->iFile, this->iBase, cgZneID, "CoordinateY",
			CG_RealDouble, &irmin, &irmax, cgZne.getYCoords());
		cg_coord_read(this->iFile, this->iBase, cgZneID, "CoordinateZ",
			CG_RealDouble, &irmin, &irmax, cgZne.getZCoords());

	}

	return true;
}

// ************** Ghost preparations ***********************


// assuming that every element-to-element connection is face-2-face connection
// number of direct ghosts is just number of elems in all connections.
// this function is needed to load cells while Ghost Manager is not initialized yet
cgsize_t t_CGNSContext::getNumOfGhostsForZone(int cgZoneID) const {

	cgsize_t nneig_cells = 0;

	const t_CGNSZone& cgZne = cgZones[cgZoneID-1];

	for (auto c : cgZne.getConns()) {
		nneig_cells+= c->get_buf().nCols;
	}

	return nneig_cells;

};




