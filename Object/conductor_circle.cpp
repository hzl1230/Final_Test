
#include "conductor_circle.h"
#include <cmath>


/* ---------------- Begin Public Methods ---------------- */

/* ------------------------------------------------------- */

bool ConductorCircle::is_inside(Real x, Real y, Real z) const
{
  // circle only for 2D
  Real dx = x-center[0], dy = y-center[1];
  bool is_inside = sqrt(dx*dx + dy*dy) <= radius;
  return is_inside;
}

/* ------------------------------------------------------- */

bool ConductorCircle::is_inside(const Vector3& point) const
{
  return is_inside(point[0], point[1], point[2]);
}

/* ------------------------------------------------------- */

bool ConductorCircle::intersect(const Vector3& pbeg,
                                const Vector3& pend,
                                std::vector<Vector3>& pintsec) const
{
  bool is_intersected = false;

  return is_intersected;
}


/* ------------------------------------------------------- */

bool ConductorCircle::intersect(const Vector3& pbeg,
                                const Vector3& pend) const
{
  std::vector<Vector3> pintsec;
  return intersect(pbeg, pend, pintsec);
}



/* ---------------- End Public Methods ---------------- */
