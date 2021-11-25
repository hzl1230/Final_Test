#ifndef _BEAM_H
#define _BEAM_H

#include "inject.h"

class BeamDef
{
  public:
    // constructor
    BeamDef(int sid, Real n, Real T, Real v[3], Real c[3], Real d[3], Real wx=1, Real wy=0)
      : specid(sid), ndens(n), temp(T), vel{v[0], v[1], v[2]},
        width_x(wx), width_y(wy), center{c[0], c[1], c[2]}, fnorm{d[0], d[1], d[2]}
        { }
    // wx = 1 by default (given in param_particle.cpp)
    // wy = 0 by default (given in param_particle.cpp)
    int specid;
    Real ndens;
    Real temp;
    Real vel[3];
    Real width_x;
    Real width_y;
    Real center[3];   // beam center
    Real fnorm[3];    // face normal vector
};

class Beam : public Inject
{
  public:
    // constructors
    Beam(int, const class BeamDef* const&, const class SpeciesDef* const&);

    void gen_particles(Real, std::vector<Particle>&);

    void print() const;

  private:
    typedef void (Beam::*GenPartFn)(Real, std::vector<Particle>&);
    GenPartFn ptr_gen_particles;

    // private methods
    void init_2d();
    void init_3d();
    void init_axi();

    void init_plane_2d_axi();
    void init_particle_cache_arr(Real);

    void gen_particles_2d(Real, std::vector<Particle>&);
    void gen_particles_3d(Real, std::vector<Particle>&);
    void gen_particles_axi(Real, std::vector<Particle>&);

    void inject_beam_particles(std::vector<Particle>&);

  protected:
    Real width_x;
    Real width_y;
    Real area;
    Real center[3];   // beam center

  private:
    int cache_size;
    int cache_cursor;
    std::vector<int> cache_count;
    std::vector<Particle> particle_cache_arr;
};

#endif
