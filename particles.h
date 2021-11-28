#ifndef _PARTICLES_H
#define _PARTICLES_H

#include <vector>
#include <random>
#include "espic_type.h"
#include "espic_math.h"

class Particle {
  public:
    // constructors
    // default constructor
    Particle() :
      pos_ {0, 0, 0},
      vel_ {0, 0, 0}, 
      nu_(0), losten(0)
      { }

    Particle(Real x, Real vx, Real y, Real vy, Real z=0., Real vz=0.) :
      pos_ {x, y, z},
      vel_ {vx, vy, vz},
      nu_(0), losten(0)
      { }
    
    Particle(Real x, Real y, Real z) :
      pos_ {x, y, z},
      vel_ {0, 0, 0},
      vel_r {0, 0, 0},
      nu_ (0), losten(0)
      { }

    // copy constructor
    Particle(const Particle& other) :
      pos_{other.pos_[0], other.pos_[1], other.pos_[2]},
      vel_{other.vel_[0], other.vel_[1], other.vel_[2]},
      vel_r{other.vel_r[0], other.vel_r[1], other.vel_r[2]},
      nu_(other.nu_), losten(other.losten)
      { }
// assignment operator
    Particle& operator=(const Particle& rhs) {
      pos_[0] = rhs.pos_[0]; pos_[1] = rhs.pos_[1]; pos_[2] = rhs.pos_[2];
      vel_[0] = rhs.vel_[0]; vel_[1] = rhs.vel_[1]; vel_[2] = rhs.vel_[2];
      vel_r[0] = rhs.vel_r[0]; vel_r[1] = rhs.vel_r[1]; vel_r[2] = rhs.vel_r[2];
      nu_ = rhs.nu_;
      losten = rhs.losten;
      return *this;
    }

    void gen_relative_vel(const Real vth)
    {
      Real vxb, vyb, vzb;
      VelBoltzDistr(vth, vxb, vyb, vzb);
      vel_r[0] = vel_[0] - vxb;
      vel_r[1] = vel_[1] - vyb;
      vel_r[2] = vel_[2] - vzb;
    }

    Real& x()  { return pos_[0]; }
    Real& y()  { return pos_[1]; }
    Real& z()  { return pos_[2]; }
    Real& vx() { return vel_[0]; }
    Real& vy() { return vel_[1]; }
    Real& vz() { return vel_[2]; }
    Real& vxr() { return vel_r[0]; }
    Real& vyr() { return vel_r[1]; }
    Real& vzr() { return vel_r[2]; }
    const Real& x()  const { return pos_[0]; }
    const Real& y()  const { return pos_[1]; }
    const Real& z()  const { return pos_[2]; }
    const Real& vx() const { return vel_[0]; }
    const Real& vy() const { return vel_[1]; }
    const Real& vz() const { return vel_[2]; }

    const Real velsqr() { return 0.5*(vel_[0]*vel_[0] 
                      + vel_[1]*vel_[1] + vel_[2]*vel_[2]); }
    const Real rel_velsqr() { return 0.5*(vel_r[0]*vel_r[0] 
                      + vel_r[1]*vel_r[1] + vel_r[2]*vel_r[2]); }
    std::vector<Real>& nu() { return nu_; }
    const std::vector<Real>& nu() const { return nu_; }
    void update_nu(const std::vector<Real>& nu) { nu_ = nu; }
    Real& lost() { return losten; }

    Real* pos() { return pos_; }
    Real* vel() { return vel_; }

  private:
    Real pos_[3];
    Real vel_[3];
    Real vel_r[3];
    std::vector<Real> nu_;
    Real losten;
    
};

// help function
inline Particle make_particle(Real x, Real y, Real vx, Real vy)
{
  return Particle(x, vx, y, vy);
}

inline Particle make_particle(Real x, Real y, Real z, Real vx, Real vy, Real vz, Real nu)
{
  return Particle(x, vx, y, vy, z, vz);
}

class Particles {
  friend class Particle;

  public:

    typedef std::size_t size_type;

    // default constructor
    Particles();

    Particles(const std::vector<Real>&, const std::vector<Real>&, const std::vector<Real>&,
              const std::vector<Real>&, const std::vector<Real>&, const std::vector<Real>&);

    // copy constructor
    Particles(const Particles&);

    // desctructor
    ~Particles();

    // public methods
    size_type size() const { return nparticles; }

//     std::vector<Real>& x()  { return pos_x; }
//     std::vector<Real>& y()  { return pos_y; }
//     std::vector<Real>& z()  { return pos_z; }
//     std::vector<Real>& vx() { return vel_x; }
//     std::vector<Real>& vy() { return vel_y; }
//     std::vector<Real>& vz() { return vel_z; }
    // returning x, y, z, vx, vy, vz as an array is used for sum_field
    const std::vector<Real>& x() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = data[ip].x();
      return scalar;
    }
    const std::vector<Real>& y() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = data[ip].y();
      return scalar;
    }
    const std::vector<Real>& z() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = data[ip].z();
      return scalar;
    }
    const std::vector<Real>& vx() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = data[ip].vx();
      return scalar;
    }
    const std::vector<Real>& vy() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = data[ip].vy();
      return scalar;
    }
    const std::vector<Real>& vz() {
      resize_scalar();
      for (size_type ip = 0; ip < size(); ip++) scalar[ip] = data[ip].vz();
      return scalar;
    }
    const std::vector<Real>& get_particles_energy() 
    {
      resize_scalar();
      for(size_type ip = 0; ip < size(); ip++) {
        Real vx2, vy2, vz2;
        vx2 = data[ip].vx() * data[ip].vx();
        vy2 = data[ip].vy() * data[ip].vy();
        vz2 = data[ip].vz() * data[ip].vz();
        scalar[ip] = 0.5 * (vx2 + vy2 + vz2);
      }
      return scalar;
    }

    // reserve space for storing n particles
    void reserve(size_type n);

    // append one particle to the end
    void append(const Particle&);
    // append a list of particles
    void append(const std::vector<Particle>&);
    // append particles given bewteen two iterators
    void append(std::vector<Particle>::const_iterator, std::vector<Particle>::const_iterator);
    void append(std::vector<Particle>::iterator, std::vector<Particle>::iterator);
    // append a Particles instance
    void append(const Particles&);
    void append(Particles*);

    void particles_shuffle();
    void get_sub_particles(size_type , Particles&);

    // "particles[id] = particle"
    // overwrite particles[id] with the new particle
    void overwrite(size_type id, const Particle& particle) {
#ifdef DEBUG
      assert(id < size());
#endif
//       pos_x[id] = particle.x();
//       pos_y[id] = particle.y();
//       pos_z[id] = particle.z();
//       vel_x[id] = particle.vx();
//       vel_y[id] = particle.vy();
//       vel_z[id] = particle.vz();
      data[id] = particle;
    }

    // particles[id1] = particles[id2]
    void overwrite(size_type id1, size_type id2) {
#ifdef DEBUG
      assert(id1 < size() && id2 < size());
#endif
//       pos_x[id1] = pos_x[id2];
//       pos_y[id1] = pos_y[id2];
//       pos_z[id1] = pos_z[id2];
//       vel_x[id1] = vel_x[id2];
//       vel_y[id1] = vel_y[id2];
//       vel_z[id1] = vel_z[id2];
      data[id1] = data[id2];
    }

    // "particle = particles[id]"
    // fetch particles[id]'s property and store in particle
    void fetch(size_type id, Particle& particle) {
#ifdef DEBUG
      assert(id < size());
#endif
//       particle.x() = pos_x[id];
//       particle.y() = pos_y[id];
//       particle.z() = pos_z[id];
//       particle.vx() = vel_x[id];
//       particle.vy() = vel_y[id];
//       particle.vz() = vel_z[id];
      particle = data[id];
    }

    Particle& operator[] (size_type i) { return data[i]; }
    const Particle& operator[] (size_type i) const { return data[i]; }
    Particle& at(size_type i) { return data[i]; }
    const Particle& at(size_type i) const { return data[i]; }

    // erase particle with id 
    void erase(size_type id);
    // erase n particles starting with id 
    void erase(size_type id, size_type n);

    // pop_back n particles from the end of array
    void pop_back(size_type n);

    // write particle info to restart file
    void write_restart(FILE*);

    void read_restart_and_store(size_type np, FILE*);

  private:
    size_type nparticles;
//     std::vector<Real> pos_x;
//     std::vector<Real> pos_y;
//     std::vector<Real> pos_z;
//     std::vector<Real> vel_x;
//     std::vector<Real> vel_y;
//     std::vector<Real> vel_z;
    std::vector<Particle> data;
    std::vector<Real> scalar;

    void resize_scalar() {
      if(scalar.size() != size()) scalar.resize(size());
    }
};

inline void RelativeVelocity(Particle& pt, Real vxb, Real vyb, Real vzb)
{
    pt.vxr() = pt.vx() - vxb;
    pt.vyr() = pt.vy() - vyb;
    pt.vzr() = pt.vy() - vzb;
    // energy = 0.5*(pt.vxr()*pt.vxr() + pt.vyr()*pt.vyr() + pt.vzr()*pt.vzr());
}

inline Real get_energy(Real vx, Real vy, Real vz) 
{
    return 0.5*(vx*vx + vy*vy + vz*vz);
}

inline Real Pcoll(Real nu, Real dt) 
{
    return 1 - exp(-nu*dt);
} 

inline Real RG01()
{
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_real_distribution<Real> RDr(0.,1.);
    return RDr(g);
}



#endif
