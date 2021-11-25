#ifndef _AMBIENT_H
#define _AMBIENT_H

#
#include "espic_type.h"

class AmbientDef {
  public:
    AmbientDef(int spid, Real n, Real T, Real v[3], Real blo[3], Real bhi[3])
      : specid (spid),
        ndens (n),
        temp (T),
        vel {v[0], v[1], v[2]},
        bound_lo {blo[0], blo[1], blo[2]},
        bound_hi {bhi[0], bhi[1], bhi[2]}
    {}

  int specid;
  Real ndens;
  Real temp;
  Real vel[3];
  Real bound_lo[3];
  Real bound_hi[3];
};

class Ambient {
  public:
    Ambient(int, const class AmbientDef* const&, const class SpeciesDef* const&);

    int species_id() const { return specid; }

    Real xmin() const { return bound_lo[0]; }
    Real ymin() const { return bound_lo[1]; }
    Real zmin() const { return bound_lo[2]; }
    Real xmax() const { return bound_hi[0]; }
    Real ymax() const { return bound_hi[1]; }
    Real zmax() const { return bound_hi[2]; }

    void gen_ambient (int [3], Real [3], Real [3], Real [3], class Particles* &);

  private:
    int specid;
    Real ndens;
    Real vth;                 // (normalized) thermal speed, sqrt(2.*temp/mass)
    Real vel[3];
    Real weight;
    Real bound_lo[3];
    Real bound_hi[3];

    typedef void (Ambient::*PtrGenAmbient)(
        int [3], Real [3], Real [3], Real [3], class Particles* &);
    PtrGenAmbient ptr_gen_ambient;
    void gen_ambient_2d (int [3], Real [3], Real [3], Real [3], class Particles* &);
    void gen_ambient_3d (int [3], Real [3], Real [3], Real [3], class Particles* &);
    void gen_ambient_axi(int [3], Real [3], Real [3], Real [3], class Particles* &);
};

#endif
