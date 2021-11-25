
#include "conductor_rectangle.h"


/* ---------------- Begin Public Methods ---------------- */

/* ------------------------------------------------------- */

bool ConductorRectangle::is_inside(Real x, Real y, Real z) const
{
  // rectangle only for 2D
  bool is_inside = (x >= blo[0] && x <= bhi[0])
                && (y >= blo[1] && y <= bhi[1]);
  return is_inside;
}

/* ------------------------------------------------------- */

bool ConductorRectangle::is_inside(const Vector3& point) const
{
  return is_inside(point[0], point[1], point[2]);
}

/* ------------------------------------------------------- */

bool ConductorRectangle::intersect(const Vector3& pbeg,
                                   const Vector3& pend,
                                   std::vector<Vector3>& pintsec) const
{
  bool is_intersected = false;

  return is_intersected;
}


/* ------------------------------------------------------- */

bool ConductorRectangle::intersect(const Vector3& pbeg,
                                   const Vector3& pend) const
{
  std::vector<Vector3> pintsec;
  return intersect(pbeg, pend, pintsec);
}



/* ---------------- End Public Methods ---------------- */
