#ifndef _CONDUCTOR_RECTANGLE_H
#define _CONDUCTOR_RECTANGLE_H

#include "conductor.h"

class ConductorRectangleDef
{
  public:
    ConductorRectangleDef()
      : id(-1), type(1), flag_collectlp(0),
        epsilon(1.), phi(0.), charge(0.)
    { }

    int id;
    int type;
    int flag_collectlp;
    Real epsilon;
    Real phi;
    Real charge;
    Vector3 center;
    Vector3 blo;
    Vector3 bhi;

    bool rf;
};

class ConductorRectangle : public Conductor
{
  public:
    // constructorsLw
    explicit ConductorRectangle(const ConductorRectangleDef& cdef)
      : Conductor(cdef.center, cdef.id, cdef.type, cdef.epsilon, cdef.phi, cdef.flag_collectlp, cdef.rf),
        blo(cdef.blo), bhi(cdef.bhi)
    { }

    // return whether or not a point is inside this object
    virtual bool is_inside(Real, Real, Real = 0.) const;
    virtual bool is_inside(const Vector3&) const;

    // find intersection points between
    // a line segment (defined by two points) and object
    virtual bool intersect(const Vector3&, const Vector3&) const;
    virtual bool intersect(const Vector3&, const Vector3&, std::vector<Vector3>&) const;

  private:
    Vector3 blo;    // boundary lo
    Vector3 bhi;    // boundary hi
};

#endif
