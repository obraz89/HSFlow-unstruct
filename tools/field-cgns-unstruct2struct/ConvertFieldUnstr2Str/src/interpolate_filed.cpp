#include "interpolate_field.h"

// get index of structured zone for the unstructured zone with index iZoneU
// detect by zone name
// TODO: correct way to detect?
int getZoneIndStruct(int iZoneU) {

	int iZoneS;

	std::string nameU = G_DomUnst.Zones[iZoneU].getName();

	for (int iZoneS = 0; iZoneS < G_Domain.nZones; iZoneS++) {

		if (nameU.compare(G_Domain.Zones[iZoneS].szName)==0)
			return iZoneS;

	}

	hsLogError("Failed to find struct zone with name %s", nameU.c_str());
	return 0;
}

// direction of snake moving in xy plane
// < (head)
// * (tail)
//<_____________ 
//______________|
//|_____________              
//*_____________|
struct t_DirXY {
	const int Nx;
	const int Ny;
	int i, j;
	bool MoveRight;
	t_DirXY() = delete;
	t_DirXY(int a_Nx, int a_Ny) : Nx(a_Nx), Ny(a_Ny), MoveRight(true) { i = 1; j = 1;}

	void move() {

		if (i == 1) {
			if (MoveRight) {
				i += 1;
				return;
			}
			else {
				MoveRight = true;
				j += 1;
				return;
			}
		}

		if (i > 1 && i < Nx) {
			i += MoveRight ? 1 : -1; 
			return;
		}

		if (i == Nx) {
			if (MoveRight) {
				MoveRight = false;
				j += 1;
				return;
			}
			else {
				i -= 1;
				return;
			}
		}

	}
};

// two domains: 
// G_DomUnst - unstructured
// G_Domain - structured
// grids must be identical (vertex2vertex)
// copy unstructured vertex values into structured field
void interpolate_match1to1() {

	const t_Dom5& G_DU = G_DomUnst;

	for (int iZoneU = G_DU.iZneMPIs; iZoneU <= G_DomUnst.iZneMPIe; iZoneU++) {

		int iZoneS = getZoneIndStruct(iZoneU);

		const t_Zone& zneU = G_DU.Zones[iZoneU];
		TZone& zneS = G_Domain.Zones[iZoneS];

		const int nVertsU = zneU.getnVerts();
		const int nVertsS = zneS.nodes_count();

		if (nVertsU != nVertsS)
			hsLogError("Vertex size mismatch: Unstruct:%d, Struct:%d", nVertsU, nVertsS);

		hsLogMessage("dims: nx=%d, ny=%d, nz=%d", zneS.nx, zneS.ny, zneS.nz);

		// "snake" iterations in plane k=const of structured block
		for (int ks = 1; ks <= zneS.nz; ks++) {
			// 1) starting vertex
			int nx_c = 1;
			int ny_c = 1;
			int nz_c = ks;

			Tcoord3D coord_c = zneS.coord<true>(nx_c, ny_c, nz_c);
			hsLogMessage("Start slice iteration, position (k=%d): x=%lf, y=%lf, z=%lf",ks, coord_c.x, coord_c.y, coord_c.z);

			t_Vec3 rc(coord_c.x, coord_c.y, coord_c.z);

			t_Vec3 dr;
			const double TOL_DR = 1.0e-09;
			int iVertC;
			for (int iVert = 0; iVert < nVertsU; iVert++) {
				dr = zneU.getVert(iVert).xyz - rc;
				if (dr.norm() < TOL_DR) {
					iVertC = iVert;
					break;
				}
			}
			hsLogMessage("iVert start:%d", iVertC);
			// move like a snake (in Nokia 3310)
			t_DirXY dir(zneS.nx, zneS.ny);
			int iVert = iVertC;
			t_PrimVars pv; 
			int IdxGlobStruct;
			const t_Vert* pVert, *pVertNeig;
			const t_Cell* pCell;
			Tcoord3D rc_coord;
			for (int iNode = 1; iNode <= zneS.nx*zneS.ny; iNode++) {

				pv = G_DU.calcVertexPV(iZoneU, iVert);
				IdxGlobStruct = zneS.globRealInd(dir.i, dir.j, ks);

				zneS.U[IdxGlobStruct + 0] = pv.getU();
				zneS.U[IdxGlobStruct + 1] = pv.getV();
				zneS.U[IdxGlobStruct + 2] = pv.getW();
				zneS.U[IdxGlobStruct + 3] = pv.getP();
				zneS.U[IdxGlobStruct + 4] = pv.getT();
				// break if this is the last point
				if (iNode == zneS.nx*zneS.ny) break;
				
				dir.move();
				rc_coord = zneS.coord<true>(dir.i, dir.j, ks);
				rc.set(rc_coord.x, rc_coord.y, rc_coord.z);

				// calculate new pair of vertices
				pVert = zneU.getpVert(iVert);
				bool found_neig = false;

				for (int iNeig = 0; iNeig < pVert->NNeigCells; iNeig++) {

					pCell = pVert->pNeigCells[iNeig];

					for (int ind = 0; ind < pCell->NVerts; ind++) {

						pVertNeig = pCell->pVerts[ind];
						dr = pVertNeig->xyz - rc;

						if (dr.norm() < TOL_DR) {
							iVert = pVertNeig->Id;
							found_neig = true;
							break;
						}
					}
				}
				if (!found_neig) {
					hsLogError("Failed to find unstruct vert:");
					hsLogError("ZoneStruct_id=%d, ind_struct=[%d, %d, %d]",
						iZoneS, dir.i, dir.j, ks);
				} 

			}
		}

	}
	// find centeral vertex

};