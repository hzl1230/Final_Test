#ifndef _COLL_
#define _COLL_
#
#include "cross_section.h"
#include "particles.h"
#include "espic_type.h"
#include "espic_math.h"

typedef std::vector<std::vector<Particle>> CollProd;
extern Real kTe0;
using namespace ESPIC;

class Collisionpair {
public:    // In class all velocity except for the Update part are relative velocity 
    Collisionpair(Particle& particle, Real m1, Real m2, Real vtb);

    ~Collisionpair();

    void ParticleElasticCollision(); 

    void ParticleExcitatinCollision(Real th);

    void ParticleIonizationCollision(Real th);

    void ParticleIsotropicCollision();

    void ParticleBackwardCollision();

    std::vector<Particle>& ion_products() { return product_arr; }

private:
    void FindEulerAngle();

    void EjectElectronReaction(Real chi_, Real eta_, Real vel_, Particle& particle);

    void UpdateParticleVelInfo();

    void EjectIonReaction(Particle& particle);


    Particle& pt;
    const Real mass;
    const Real vth;
    Real gx, gy, gz, gyz, g;    // relative-velocity
    const Real F1, F2;
    
    Real chi, eta;
    Real wx, wy, wz;
    Real st, ct, cp, sp, stcp, stsp;
    std::vector<Particle> product_arr;
    Real energy;
    // std::vector<Real> velbuffer;

};
#endif