#ifndef _COLL_
#define _COLL_
#
#include "particles.h"
#include "espic_type.h"
#include "espic_math.h"
#include "species.h"

typedef std::vector<std::vector<Particle>> CollProd;
extern Real kTe0;
using namespace ESPIC;

class Collisionpair {
public:    // In class all velocity except for the Update part are relative velocity 
    Collisionpair(Real m1, Real m2, Real vtb);

    ~Collisionpair();

    void gen_collision_info(const Particle& ip, const Real* vr, const Real vel, const Real energy, const bool e_col_flag);

    void ElasticScatter(Particle& particle); 

    void ExcitatinCollision(const Real th, Particle& particle);

    void IonizationCollision(const Real th, Particle& particle,
                             const std::vector<int>& prodid, 
                             std::vector<Species*>& sp_arr);

    void DissociationAttechment(const Real th, 
                                Particle& particle,
                                const std::vector<int>& prodid, 
                                std::vector<Species*>& sp_arr);
    
    void DetachmentCollision(const Real th, 
                    Particle& particle,
                    const std::vector<int>& prodid,
                    std::vector<Species*>& sp_arr);

    void IsotropicCollision(Particle& particle);

    void BackwardCollision(Particle& particle);

    void ChargeExchange(Particle& particle);

    void CoulombScatter(Particle& ip, Particle& tp);

    void ParticleRecombination(Particle& ip, Particle& tp);

    void ElectronImpactDetechment(const Real th, Particle& ip,
                                  const std::vector<int>& prodid, 
                                  std::vector<Species*>& sp_arr);

    

    // std::vector<Particle>& ion_products() { return product_arr; }
    const Real vx() { return v_new[0]; }
    const Real vy() { return v_new[1]; }
    const Real vz() { return v_new[2]; }

private:
    void FindEulerAngle();

    void EjectElectronReaction(Real chi_, Real eta_, Real vel_, 
                               Real& vx, Real& vy, Real& vz);

    void UpdateParticleVelInfo(Particle& particle);

    // void EjectIonReaction(Real& vx, Real& vy, Real& vz);

    const Real mr;
    const Real vth;
    const Real F1, F2;
    // Particle& particle;
    Real gx, gy, gz, gyz;    // relative-velocity
    Real g;
    Real energy;

    Real chi, eta;
    Real wx, wy, wz;
    Real st, ct, cp, sp, stcp, stsp;
    Real v_new[3];
};
#endif