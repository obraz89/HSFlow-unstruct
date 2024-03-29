#include "flux_ns.h"

#include "flow_model_perfect_gas.h"

#include "flow_params.h"


//******************************************Flux
// all input vectors must be in global reference frame

void calcNSViscFlux(const t_Vec3& Norm, const t_PrimVars& PV, const t_Mat<3, 5> GradRUVWT, t_VecConsVars& flux) {

	const double
		vx = PV.getU(), vy = PV.getV(), vz = PV.getW(), p = PV.getP(), T = PV.calcT();

	const double
		vx_dx = GradRUVWT[0][1], vx_dy = GradRUVWT[1][1], vx_dz = GradRUVWT[2][1];

	const double
		vy_dx = GradRUVWT[0][2], vy_dy = GradRUVWT[1][2], vy_dz = GradRUVWT[2][2];

	const double
		vz_dx = GradRUVWT[0][3], vz_dy = GradRUVWT[1][3], vz_dz = GradRUVWT[2][3];

	const double
		T_dx = GradRUVWT[0][4], T_dy = GradRUVWT[1][4], T_dz = GradRUVWT[2][4];


	const double div = vx_dx + vy_dy + vz_dz;

	const double s_xx = 2 * vx_dx - 2. / 3 * div; 
	const double s_xy = vx_dy + vy_dx;      
	const double s_xz = vx_dz + vz_dx;      
	const double s_yy = 2 * vy_dy - 2. / 3 * div; 
	const double s_yz = vy_dz + vz_dy;      
	const double s_zz = 2 * vz_dz - 2. / 3 * div; 

	double mju_Re = calcViscosity(T)/G_FreeStreamParams.getRe();
	double lam_Re = calcThermalConductivity(T)/G_FreeStreamParams.getRe();

	t_VecConsVars E, G, F;
	
	// flux x
	E[0] = 0.0;
	E[1] = -mju_Re*s_xx;
	E[2] = -mju_Re*s_xy;
	E[3] = -mju_Re*s_xz;
	E[4] = -mju_Re*(vx*s_xx + vy*s_xy + vz*s_xz) - lam_Re*T_dx;

	// flux y
	G[0] = 0.0;
	G[1] = -mju_Re*s_xy;
	G[2] = -mju_Re*s_yy;
	G[3] = -mju_Re*s_yz;
	G[4] = -mju_Re*(vx*s_xy + vy*s_yy + vz*s_yz) - lam_Re*T_dy;

	// flux z
	F[0] = 0.0;
	F[1] = -mju_Re*s_xz;
	F[2] = -mju_Re*s_yz;
	F[3] = -mju_Re*s_zz;
	F[4] = -mju_Re*(vx*s_xz + vy*s_yz + vz*s_zz) - lam_Re*T_dz;

	for (int i = 0; i < NConsVars; i++) {
		flux[i] = Norm[0] * E[i] + Norm[1] * G[i] + Norm[2] * F[i];
	}

}



