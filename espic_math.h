#ifndef ESPIC_MATH_H
#define ESPIC_MATH_H

#include <functional>
#include <limits>
#include <random>
#include <cmath>

#include "espic_type.h"
namespace ESPIC{
  // constants
  const Real PI  = 2.0*asin(1.0);
  const Real PI2 = 4.0*asin(1.0);
  const Real SMALLREAL = std::numeric_limits<Real>::min();

  // Newton-Raphson method
  const int maxit_nr = 100;
  const Real rtol_nr = 1e-12;
  Real newton_raphson(const std::function< Real(Real) >&, 
                      const std::function< Real(Real) >&);

  static std::vector<Real> velbuffer; // Store Random Vel bd

  /* class of random number generators */
  /* generate uniform distributions in [0, 1)*/
  class Random {
    public :
      /* Constructor */
      Random() : p{0.0, 1.0} 
      {
        if (seed < 0) {
          auto rs = std::bind(std::uniform_int_distribution<> (), 
                              std::default_random_engine ());
          seed = rs();
        }
        else {
          auto rs = std::bind(std::uniform_int_distribution<> (), 
                              std::default_random_engine (seed));
          seed = rs();
        }
        
        rf = std::bind(std::uniform_real_distribution<Real> {p},
            std::default_random_engine(seed));
      }

      Real operator() () const { return rf(); }     // generator a random number

      Real uniform_dist() const { return rf(); }

      Real normal_dist_factor() const {
        Real r = rf();
        while(r <= SMALLREAL) {
          r = rf();
        }
        return sqrt(-log(r));
      }

      Bigint get_seed() const   { return seed; }

      static void set_seed(Bigint s) { seed = s; }

      

    private :
      static Bigint seed;
      std::uniform_real_distribution<Real>::param_type p;
      std::function<Real()> rf;
  };

inline void cross_prod(const Real vecA[3], const Real vecB[3], Real vecC[3])                                      
{
  double temp[3];

  temp[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
  temp[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
  temp[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];

  for(int a = 0; a < 3; a++) vecC[a] = temp[a];
}

};

extern ESPIC::Random ranf;
ESPIC::Random ranf;

inline void VelBoltzDistr(Real vth, Real& vx, Real& vy, Real& vz)
{
    
    Real vrf = vth * ranf.normal_dist_factor();
    Real trf = ESPIC::PI2 * ranf();
    vx = vrf * cos(trf);
    vy = vrf * sin(trf);
    
    vrf = vth * ranf.normal_dist_factor();
    trf = ESPIC::PI2 * ranf();
    vz = vrf * cos(trf);
    // velbuffer.emplace_back(vrf*sin(trf));
}

#endif
