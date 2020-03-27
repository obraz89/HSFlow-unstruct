#include "matrix_small.h"

t_Vec3 operator*(double val, const t_Vec3& vec) {
	t_Vec3 ret;
	(t_Vec<3>)ret= vec*val;
	return ret;
};