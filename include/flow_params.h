#pragma once

#include "PluginBase.h"

#include "matrix_small.h"

#include <cmath>

#include "settings.h"

#define STR_FLOW_PARAMS "flow_params"

class t_FlowParamsFreeStream : public TPlugin {
protected:
	std::string name;

	double Mach;

	double TinfDim;

	// aoa in radians
	double AoARad;

	// aos in radians
	double AoSRad;

	// unit vector of velocity direction at infinity
	t_Vec3 UInf;

	// non-dimensional Reynolds number via values at infinity & 
	// reference Length used to non-dim grid 
	double Re;
public:
	t_FlowParamsFreeStream() :name(STR_FLOW_PARAMS) { default_settings(); }
	// implement TPugin
	std::string get_description() const { return std::string("flow params"); };
	std::string get_name() const { return name; };
	void default_settings() {
		TPluginParamsGroup g("", "gas-dynamic values at infinity");
		g.add("Mach", 2.0, "Mach number");
		g.add("Tinf", 300.0, "Temperature");
		g.add("AoA", 0.0, "Angle of attack (degrees)");
		g.add("AoS", 0.0, "Angle of glide (Skol'zhenie) (degrees)");
		g.add("Re", 1.0e+06, "Reynolds number");
		_mapParamsGrps.emplace(g.get_name(), g);
	};

	void init(std::string& ini_data, const std::string& spec) {

		TPlugin::init(ini_data, spec);

		const TPluginParamsGroup& g = get_settings_grp("");

		// copied PI from std::M_PI
		double PI = 3.14159265358979323846;

		Mach = g.get_real_param("Mach");

		TinfDim = g.get_real_param("Tinf");

		double aoa_deg = g.get_real_param("AoA");

		AoARad = PI * aoa_deg / 180.0;

		double aos_deg = g.get_real_param("AoS");

		AoSRad = PI * aos_deg / 180.0;

		double UxDir = cos(AoARad) * cos(AoSRad);

		double UyDir = sin(AoARad);

		double UzDir = cos(AoARad) * sin(AoSRad);
		// free stream velocity vector
		switch (g_genOpts.nonDimType.Val)
		{
		case t_NonDimType::FreeStreamVelo:
			UInf.set(UxDir, UyDir, UzDir);
			break;
		case t_NonDimType::FreeStreamSoundSpeed:
			UInf.set(Mach * UxDir, Mach * UyDir, Mach * UzDir);
			break;
		default:
			hsLogError("Flow Params: unknown nondim type");
			break;
		}


		Re = g.get_real_param("Re");
	};

	const t_Vec3& getUInf() const { return UInf; }
	double getMach() const { return Mach; }
	double getTinfDim() const { return TinfDim; }
	double getRe() const { return Re; }

};

extern t_FlowParamsFreeStream G_FreeStreamParams;

//-----------------------------------------------------------------------------
