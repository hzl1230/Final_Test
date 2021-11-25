#ifndef _OBJECT_H
#define _OBJECT_H

#include "../espic_type.h"

class Object {
  public:
    // constructors
    explicit Object(int _id = -1)
      : center{0., 0., 0.} {
        if (_id > 0) current_id = id = _id;
        else id = get_next_id();
      }
      
    Object(Real c[3], int _id = -1)
      : center{c[0], c[1], c[2]} {
        if (_id > 0) current_id = id = _id;
        else id = get_next_id();
    }

    Object(const Vector3& c, int _id = -1)
      : center{c[0], c[1], c[2]} {
        if (_id > 0) current_id = id = _id;
        else id = get_next_id();
    }

    // copy constrcutor (next_id unchanged)
    Object(const Object& orig)
      : id(orig.id),
        center{orig.center[0], orig.center[1], orig.center[2]}
    { }

    // destructor
    virtual ~Object() { }

    // public methods
    static int get_next_id() { return ++current_id; }

    int get_id() const { return id; }

    // return whether or not a point is inside this object
    virtual bool is_inside(Real, Real, Real = 0.) const = 0;
    virtual bool is_inside(const Vector3&) const = 0;

    // find intersection points between
    // a line segment (defined by two points) and object
    virtual bool intersect(const Vector3&, const Vector3&) const = 0;
    virtual bool intersect(const Vector3&, const Vector3&,
                           std::vector<Vector3>&) const = 0;

  protected:
    static int current_id;
    int id;
    Real center[3];
};

#endif
