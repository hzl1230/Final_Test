#include "collision.h"
#include "espic_type.h"
#include <fstream>

// In class all velocity except for the Update part are relative velocity 
Collisionpair::Collisionpair(Real m1, Real m2, Real vtb)
: mr(m1*m2/(m1+m2)), 
vth(vtb),
F1(m1/(m1+m2)), 
F2(m2/(m1+m2))
{ }

Collisionpair::~Collisionpair()
{ }

void Collisionpair::gen_collision_info(const Particle& particle, const Real* vr, const Real vel, 
                                       const Real en, const bool e_col_flag)
{
    energy = en; g = vel;
    gx = vr[0]; gy = vr[1]; gz = vr[2];
    if (!e_col_flag) {
        Real vxb, vyb, vzb;
        vxb = particle.vx() - gx;
        vyb = particle.vy() - gy;
        vzb = particle.vz() - gz;
        wx = F1 * particle.vx() + F2 * vxb;
        wy = F1 * particle.vy() + F2 * vyb;
        wz = F1 * particle.vz() + F2 * vzb;
    } else {
        wx = F1 * vr[0];
        wy = F1 * vr[1];
        wz = F1 * vr[2];
    }
}

void Collisionpair::ElasticScatter(Particle& particle) 
{ 
    chi = acos(1.0 - 2.0*ranf());
    eta = ESPIC::PI2 * ranf();
    Real sc(sin(chi)), cc(cos(chi));
    Real se(sin(eta)), ce(cos(eta));
    cp = gy / gyz;
    sp = gz / gyz;

    Real vx, vy, vz;
    vx = gx * cc - gyz * sc * ce;
    vy = gy * cc + gx * cp * sc * ce - g * sp * sc * se;
    vz = gz * cc + gx * sp * sc * ce + g * cp * sc * se;
    particle.vx() = wx + F2 * vx;
    particle.vy() = wy + F2 * vy;
    particle.vz() = wz + F2 * vz;
}

void Collisionpair::ExcitatinCollision(const Real th, Particle& particle) 
{
    FindEulerAngle();
    energy = fabs(energy - th);
    g = sqrt(2.0 * energy / mr);
    chi = acos(1.0 - 2.0 * ranf());
    eta = ESPIC::PI2 * ranf();
    UpdateParticleVelInfo(particle);
}


void Collisionpair::IonizationCollision(const Real th, 
                                        Particle& particle,
                                        const std::vector<int>& prodid, 
                                        std::vector<Species*>& sp_arr)
{
    Real en_ej, en_sc;
    Real g_ej, chi_ej, eta_ej;
    Real w = 10.3 / kTe0;
    energy = energy - th;
    en_ej = w * tan(ranf() * atan(0.5*energy/w));
    en_sc = fabs(energy - en_ej);
    g = sqrt(2.0 * en_sc/mr);
    g_ej = sqrt(2.0 * en_ej/mr);
    chi = acos(sqrt(en_sc / energy));
    chi_ej = acos(sqrt(en_ej / energy));
    eta = ESPIC::PI2 * ranf();
    eta_ej = eta + ESPIC::PI;
    FindEulerAngle();
    UpdateParticleVelInfo(particle);

    Real vx, vy, vz;
    Real x(particle.x()), y(particle.y()), z(particle.z());
    for(const int spid: prodid) {
        if("e" == sp_arr[spid]->name)
            EjectElectronReaction(chi_ej, eta_ej, g_ej, vx, vy, vz);
        else { VelBoltzDistr(vth, vx, vy, vz); }
    
        sp_arr[spid]->particles->append(Particle(x, y, z, vx, vy, vz));
    }
}

void Collisionpair::DetachmentCollision(const Real th, 
                               Particle& particle,
                               const std::vector<int>& prodid,
                               std::vector<Species*>& sp_arr)
{
    energy -= th;
    Real en_ej(th), g_ej, chi_ej, eta_ej;
    g_ej = sqrt(2.0 * en_ej/mr);
    chi_ej = acos(sqrt(en_ej / energy));
    eta_ej = ESPIC::PI2 + ESPIC::PI;
    
    Real vx, vy, vz;
    Real x(particle.x()), y(particle.y()), z(particle.z());
    for(const int spid: prodid) {
        if("e" == sp_arr[spid]->name)
            EjectElectronReaction(chi_ej, eta_ej, g_ej, vx, vy, vz);
        else { VelBoltzDistr(vth, vx, vy, vz); }
    
        sp_arr[spid]->particles->append(Particle(x, y, z, vx, vy, vz));
    }
}

void Collisionpair::IsotropicCollision(Particle& particle)
{
    chi = acos(1.0 - 2.0*ranf());
    eta = PI2 * ranf();
    FindEulerAngle();
    UpdateParticleVelInfo(particle);
}

void Collisionpair::BackwardCollision(Particle& particle)
{
    chi = PI;
    eta = PI2 * ranf();
    FindEulerAngle();
    UpdateParticleVelInfo(particle);
}

void Collisionpair::ChargeExchange(Particle& particle)
{
    particle.vx() = 0.0;
    particle.vy() = 0.0;
    particle.vz() = 0.0;
}

void Collisionpair::CoulombScatter(Particle& ip, Particle& tp)
{
}

void Collisionpair::DissociationAttechment(const Real th, 
                                           Particle& particle,
                                           const std::vector<int>& prodid, 
                                           std::vector<Species*>& sp_arr)
{ 
}

void Collisionpair::ElectronImpactDetechment(const Real th, Particle& ip,
                                             const std::vector<int>& prodid, 
                                             std::vector<Species*>& sp_arr)
{
    Real en_ej, en_sc;
    Real g_ej, chi_ej, eta_ej;
    Real w = 10.3 / kTe0;
    energy = energy - th;
    en_ej = w * tan(ranf() * atan(0.5*energy/w));
    en_sc = fabs(energy - en_ej);
    g = sqrt(2.0 * en_sc/mr);
    g_ej = sqrt(2.0 * en_ej/mr);
    chi = acos(sqrt(en_sc / energy));
    chi_ej = acos(sqrt(en_ej / energy));
    eta = ESPIC::PI2 * ranf();
    eta_ej = eta + ESPIC::PI;
    FindEulerAngle();
    UpdateParticleVelInfo(ip);

    Real vx, vy, vz;
    Real x(ip.x()), y(ip.y()), z(ip.z());
    for(const int spid: prodid) {
        if("e" == sp_arr[spid]->name)
            EjectElectronReaction(chi_ej, eta_ej, g_ej, vx, vy, vz);
        else { VelBoltzDistr(vth, vx, vy, vz); }
    
        sp_arr[spid]->particles->append(Particle(x, y, z, vx, vy, vz));
    }
}

void Collisionpair::FindEulerAngle()
{
    st = gyz / g;
    ct = gx / g;
    stcp = gy / g;
    stsp = gz / g;
    if (gyz == 0.) {
        sp = 0.;
        cp = 0.;
    } else {
        sp = gz / gyz;
        cp = gy / gyz;
    }
}


void Collisionpair::EjectElectronReaction(Real chi_, Real eta_, Real vel_, 
                                          Real& vx, Real& vy, Real& vz)
{
    Real sc(sin(chi_)), cc(cos(chi_));
    Real se(sin(eta_)), ce(cos(eta_));
    // Real vx, vy, vz;
    gx = vel_ * (ct * cc - st * sc * ce);
    gy = vel_ * (stcp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = vel_ * (stsp * cc + ct * sp * sc * ce + cp * sc * se);
    vx = wx + F2*gx;
    vy = wy + F2*gy;
    vz = wz + F2*gz;
}

void Collisionpair::UpdateParticleVelInfo(Particle& particle)
{
    Real sc(sin(chi)), cc(cos(chi));
    Real se(sin(eta)), ce(cos(eta));
    gx = g * (ct * cc - st * sc * ce);
    gy = g * (stcp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (stsp * cc + ct * sp * sc * ce + cp * sc * se);
    particle.vx() = wx + F2 * gx;
    particle.vy() = wy + F2 * gy;
    particle.vz() = wz + F2 * gz;
}
