#ifndef _CONDUCTOR_H
#define _CONDUCTOR_H

#include <vector>
#include <utility>
#include "object.h"

class Conductor : public Object {
  public:
    // constructors
    explicit Conductor(int _id = -1, int _type = 1, Real _epsilon = 1., Real _phi = 0., int lparticles = 0, bool is_rf = false)
      : Object(_id), type(_type), flag_collectlp(lparticles), epsilon(_epsilon), phi(_phi), charge(0.), rf(is_rf)
    { }

    Conductor(Real c[3], int _id = -1, int _type = 1, Real _epsilon = 1., Real _phi = 0., int lparticles = 0, bool is_rf = false)
      : Object(c, _id), type(_type), flag_collectlp(lparticles), epsilon(_epsilon), phi(_phi), charge(0.), rf(is_rf)
    { }
      
    Conductor(const Vector3& c, int _id = -1, int _type = 1, Real _epsilon = 1., Real _phi = 0., int lparticles = 0, bool is_rf = false)
      : Object(c, _id), type(_type), flag_collectlp(lparticles), epsilon(_epsilon), phi(_phi), charge(0.), rf(is_rf)
    { }
      
    // copy constructors
    Conductor(const Conductor& orig)
      : Object(orig), type(orig.type), flag_collectlp(orig.flag_collectlp),
        epsilon(orig.epsilon), phi(orig.phi), charge(orig.charge),
        lostparticles(orig.lostparticles)
    { }

    // desctructor
    virtual ~Conductor() { }

    // public methods
    bool is_real() const { return type > 0; }

    bool is_fixed_potential() const { return type < 2; }

    bool is_lost_particles_statistics_on() const { return flag_collectlp == 1; }

    void lostparticles_increment(int, Real);

    Real get_epsilon() const { return epsilon; }

    Real get_potential() const { return phi; }

    Real get_charge() const { return charge; }

    bool is_rf() const { return rf; }
    
    // return whether or not a point is inside this object
    virtual bool is_inside(Real, Real, Real = 0.) const = 0;
    virtual bool is_inside(const Vector3&) const = 0;

    // find intersection points between
    // a line segment (defined by two points) and object
    virtual bool intersect(const Vector3&, const Vector3&) const = 0;
    virtual bool intersect(const Vector3&, const Vector3&,
                           std::vector<Vector3>&) const = 0;

    bool scrape_particle(int, Real, const Vector3&, const Vector3&);

  protected:
    class LostParticleInfo {
      public:
        LostParticleInfo(int nspecies, Real t)
          : time(t), num_lostparticles(nspecies)
        { }

      bool is_used() const {
        bool is_used = false;
        for (auto it = num_lostparticles.cbegin(); !is_used && it != num_lostparticles.cend(); ++it)
          is_used = *it > 0;
        return is_used;
      }

      void add_one_lost(int spec_id) { ++num_lostparticles[spec_id]; }

      Real time;
      std::vector<Bigint> num_lostparticles;
    };

    int type;               // 0 - virtual, (fixed phi)
                            // 1 - real, fixed phi
                            // 2 - real, floating phi
    int flag_collectlp;     // 0 - does not count
                            // 1 - count
    Real epsilon;           // relative permittivity
    Real phi;               // potential
    Real charge;
    std::vector<LostParticleInfo> lostparticles;
    bool rf;                // true if rf is rf_source

    // private methods
    void accumulate_charge(std::vector<class Species*>&);
};

#endif
