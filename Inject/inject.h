#ifndef _INJECT_H
#define _INJECT_H

#include <vector>
#include "../particles.h"
#include "../espic_type.h"

class InjectDef {
  public:
    // constructors
    InjectDef();

    // destructor
    ~InjectDef();

    // public methods
    int num_beams() const { return static_cast<int> (beamdef_arr.size()); }

    void append(class BeamDef* b);

//     void append(class FlowDef* f);

    // members
    std::vector<class BeamDef*> beamdef_arr;
//     std::vector<class FlowDef*> flowdef_arr;

  private:

    // forbid copy
    InjectDef(const InjectDef&);
};

class Inject {  // base for particle injection
  public:
    // constructors
    Inject()
      : specid(0), ndens(0), temp(0), vel {0, 0, 0} { }

    Inject(const int sid, const Real n, const Real T, const Real v[3], Real nr=0.)
      : specid(sid), ndens(n), temp(T), vel {v[0], v[1], v[2]} , nres(nr) { }

    // desctructor
    virtual ~Inject() { }

    virtual void gen_particles(Real, std::vector<Particle>&) = 0;

    virtual void print() const = 0;

    int species() { return specid; }

  protected:
    int specid;         // species id
    Real ndens;         // number density
    Real temp;          // temperature
    Real vel[3];        // drifting velocity
    Real vth;           // thermal speed
    Real weight;
    Real sn;            // speed ratio normal to beam plane
    Real pn;            // precomputed parameter
    Real nflowrate;     // particle number flow rate
    Real nres;          // residual of injected particle
    Real tmat[3][3];    // transformation matrix

    // preform some pre-computation for part injection
    void precomputed(Real);

    // calculate number flux density
    Real calc_nflux();
    
    // generate velocity for one part to be injected
    void gen_one_vel(Real v[3]);

};

#endif
