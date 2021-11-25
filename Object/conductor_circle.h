#ifndef _CONDUCTOR_CIRCLE_H
#define _CONDUCTOR_CIRCLE_H

#include "conductor.h"

class ConductorCircleDef
{
  public:
    ConductorCircleDef()
      : id(-1), type(1), flag_collectlp(0),
        epsilon(1.), phi(0.), charge(0.), radius(0.)
    { }

    int id;
    int type;
    int flag_collectlp;
    Real epsilon;
    Real phi;
    Real charge;
    Real radius;
    Vector3 center;
};

class ConductorCircle : public Conductor
{
  public:
    // constructorsLw
    explicit ConductorCircle(const ConductorCircleDef& cdef)
      : Conductor(cdef.center, cdef.id, cdef.type, cdef.epsilon, cdef.phi, cdef.flag_collectlp),
        radius(cdef.radius)
    { }

    // return whether or not a point is inside this object
    virtual bool is_inside(Real, Real, Real = 0.) const;
    virtual bool is_inside(const Vector3&) const;

    // find intersection points between
    // a line segment (defined by two points) and object
    virtual bool intersect(const Vector3&, const Vector3&) const;
    virtual bool intersect(const Vector3&, const Vector3&, std::vector<Vector3>&) const;

  private:
    Real radius;
};

#endif
