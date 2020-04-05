#include "common_data.h"

#include "common_procs.h"

// stat(), getcwd() functions
#if defined(_WINDOWS)
#include <direct.h>  // mkdir
#define stat  _stat64
#define getcwd  _getcwd
#else
#include <unistd.h>
#endif

void ComputeTriangleAreaNormal(const t_Vec3(&pnts)[3], t_Vec3& norm, double& area) {

	t_Vec3 v1 = pnts[1] - pnts[0];
	t_Vec3 v2 = pnts[2] - pnts[0];

	norm = v1.cross(v2);
	area = 0.5*norm.norm();
	norm.normalize();

	//debug
	//hsLogMessage("Tria:Normal:(%lf, %lf, %lf), Area:%lf", norm.x, norm.y, norm.z, area);

};

// TODO: this can be computed in different ways
// here use decomposition into triangles
// if points are v1,v2,v3,v4 and center is v5
// compute via triangles (1,2,5), (2,3,5), (3,4,5), (4,1,5)
void ComputeQuadAreaNormal_v1(const t_Vec3(&pnts)[4], t_Vec3& norm, double& area) {

	t_Vec3 n_v[4];
	double s_v[4];

	t_Vec3 c;

	for (int i = 0; i < 4; i++) c += pnts[i];
	c *= 0.25;

	ComputeTriangleAreaNormal({pnts[0], pnts[1], c}, n_v[0], s_v[0]);
	ComputeTriangleAreaNormal({pnts[1], pnts[2], c}, n_v[1], s_v[1]);
	ComputeTriangleAreaNormal({pnts[2], pnts[3], c}, n_v[2], s_v[2]);
	ComputeTriangleAreaNormal({pnts[3], pnts[0], c}, n_v[3], s_v[3]);

	area = 0; for (int i = 0; i < 4; i++) area += s_v[i];
	norm.set(0, 0, 0);
	double inv_area = 1.0 / area;
	for (int i = 0; i < 4; i++) {
		double coef = s_v[i] * inv_area;
		norm += coef*n_v[i];
	} ;

	//debug
	//hsLogMessage("Quad[v1]:Normal:(%lf, %lf, %lf), Area:%lf", norm.x, norm.y, norm.z, area);

}
// here use simple quad approach
void ComputeQuadAreaNormal_v2(const t_Vec3(&pnts)[4], t_Vec3& norm, double& area) {

	t_Vec3 v1,v2;

	v1 = 0.5*(pnts[1] - pnts[0] + pnts[2] - pnts[3]);
	v2 = 0.5*(pnts[3] - pnts[0] + pnts[2] - pnts[1]);

	t_Vec3 p = v1.cross(v2);
	area = p.norm();
	norm = p; norm.normalize();

	//debug
	//hsLogMessage("Quad[v2]:Normal:(%lf, %lf, %lf), Area:%lf", norm.x, norm.y, norm.z, area);

}

void ComputeQuadAreaNormal(const t_Vec3(&pnts)[4], t_Vec3& norm, double& area) {
	ComputeQuadAreaNormal_v2(pnts, norm, area);
}

/**
 * Checks for a file existence
 *
 * @param fn file name
 *
 * @return true if exist, false otherwise
 */
bool hs_file_exists(const std::string& fn)
{
	struct stat buf;
	if (stat(fn.c_str(), &buf) == 0)
		return bool(buf.st_mode & S_IFREG);

	return false;
}
//-----------------------------------------------------------------------------

/**
 * Checks for a directory existence
 *
 * @param dn directory name
 *
 * @return true if exist, false otherwise
 */
bool hs_dir_exists(const std::string& dn)
{
	struct stat buf;
	if (stat(dn.c_str(), &buf) == 0)
		return bool(buf.st_mode & S_IFDIR);

	return false;
}
//-----------------------------------------------------------------------------


/**
 * Get file/directory modification time
 *
 * @param fn file name
 *
 * @return modification time, -1 on error
 */
time_t hs_file_mtime(const std::string& fn)
{
	struct stat buf;
	if (stat(fn.c_str(), &buf) == 0)
		return buf.st_mtime;

	return time_t(-1);
}
//-----------------------------------------------------------------------------

/**
 * Create a directory, only a last component of dirname without parents(!)
 *
 * @param dn directory name
 *
 * @return true if created, false otherwise
 */
bool hs_dir_create(const std::string& dn)
{
	if (hs_dir_exists(dn))
		return true;

	int rc = 0;

#if defined(_WINDOWS)
	rc = _mkdir(dn.c_str());
#else
	mode_t nMode = 0777; // UNIX style permissions, will be reduced by umask
	rc = mkdir(dn.c_str(), nMode);
#endif

	return rc == 0;
}
//-----------------------------------------------------------------------------



