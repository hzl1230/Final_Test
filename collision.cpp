#include "collision.h"
#include <fstream>

// In class all velocity except for the Update part are relative velocity 
Collisionpair::Collisionpair(Particle& particle, Real m1, Real m2, Real vtb)
: pt(particle), 
mass(m1*m2/(m1+m2)), vth(vtb),
gx(pt.vxr()), gy(pt.vyr()), gz(pt.vzr()),
g(sqrt(2.0 * particle.rel_velsqr())),
F1(m1/(m1+m2)), F2(m2/(m1+m2))
{
    gyz = sqrt(gy*gy + gz*gz);
    g1 = g;

    Real vxb, vyb, vzb;
    vxb = pt.vx() - gx;
    vyb = pt.vy() - gy;
    vzb = pt.vz() - gz;

    wx = F1 * pt.vx() + F2 * vxb;
    wy = F1 * pt.vy() + F2 * vyb;
    wz = F1 * pt.vz() + F2 * vzb;

    energy = pt.rel_velsqr() * mass;
}

Collisionpair::~Collisionpair()
{ }

void Collisionpair::ParticleElasticCollision() 
{ 
    std::ofstream of("ela.dat", std::ofstream::app);
    of << "before: " << pt.velsqr()*mass;
    of << "en: " << pt.rel_velsqr()*mass;

    chi = acos(1.0 - 2.0*RG01());
    eta = ESPIC::PI2 * RG01();
    Real sc(sin(chi)), cc(cos(chi));
    Real se(sin(eta)), ce(cos(eta));
    cp = gy / gyz;
    sp = gz / gyz;

    Real vx, vy, vz;
    vx = gx * cc - gyz * sc * ce;
    vy = gy * cc + gx * cp * sc * ce - g * sp * sc * se;
    vz = gz * cc + gx * sp * sc * ce + g * cp * sc * se;
    pt.vx() = wx + F2 * vx;
    pt.vy() = wy + F2 * vy;
    pt.vz() = wz + F2 * vz;

    of << " after: " << pt.velsqr()*mass << std::endl;
    of << " en: " << pt.rel_velsqr()*mass << std::endl;
}

void Collisionpair::ParticleExcitatinCollision(Real th) 
{
    std::ofstream of("exc.dat", std::ofstream::app);
    of << "th: " << th;
    of << " before: " << pt.velsqr()*mass;
    FindEulerAngle();
    energy = fabs(energy - th);
    g1 = sqrt(2.0 * energy / mass);
    chi = acos(1.0 - 2.0 * RG01());
    eta = ESPIC::PI2 * RG01();
    UpdateParticleVelInfo();
    pt.lost() += th;

    of << " after: " << pt.velsqr()*mass << std::endl;
    of << std::endl;
    of.close();
}

void Collisionpair::ParticleIonizationCollision(Real th)
{
    Real en_ej, en_sc;
    Real g_ej, chi_ej, eta_ej;
    Real w = 10.3 / kTe0;

    std::ofstream of("ion.dat", std::ofstream::app);

    energy = fabs(energy - th);
    en_ej = w * tan(RG01() * atan(0.5*energy/w));
    en_sc = fabs(energy - en_ej);
    g1 = sqrt(2.0 * en_sc/mass);
    g_ej = sqrt(2.0 * en_ej/mass);
    chi = acos(sqrt(en_sc / energy));
    chi_ej = acos(sqrt(en_ej / energy));
    eta = ESPIC::PI2 * RG01();
    eta_ej = eta + ESPIC::PI;

    of << "th: " << th ;
    of << " before: " << pt.velsqr()*mass;

    Particle e_ej = Particle(pt.x(), pt.y(), pt.z());
    Particle p_ej = Particle(pt.x(), pt.y(), pt.z());

    FindEulerAngle();
    UpdateParticleVelInfo();

    of << " after: " << pt.velsqr()*mass << std::endl;

    EjectElectronReaction(chi_ej, eta_ej, g_ej, e_ej);
    product_arr.emplace_back(std::move(e_ej));
    EjectIonReaction(p_ej);
    product_arr.emplace_back(std::move(p_ej));

    of << "Eng sc: " << en_sc << " "
       << "Eng ej: " << en_ej << std::endl;
    of << "F1: " << F1 << " "
       << "F2: " << F2 << std::endl;
    of << "sc: " << pt.velsqr()*mass << " "
       << "ej: " << e_ej.velsqr()*mass;
       
       
    pt.lost() += th;
    of << std::endl;
    of.close();
}

void Collisionpair::ParticleIsotropicCollision()
{
    chi = acos(1.0 - 2.0*RG01());
    eta = PI2 * RG01();
    FindEulerAngle();
    UpdateParticleVelInfo();
}

void Collisionpair::ParticleBackwardCollision()
{
    chi = PI;
    eta = PI2 * RG01();
    FindEulerAngle();
    UpdateParticleVelInfo();
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


void Collisionpair::EjectElectronReaction(Real chi_, Real eta_, Real vel_, Particle& particle)
{
    Real sc(sin(chi_)), cc(cos(chi_));
    Real se(sin(eta_)), ce(cos(eta_));
    // Real vx, vy, vz;
    gx = vel_ * (ct * cc - st * sc * ce);
    gy = vel_ * (stcp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = vel_ * (stsp * cc + ct * sp * sc * ce + cp * sc * se);
    particle.vx() = wx + F2*gx;
    particle.vy() = wy + F2*gy;
    particle.vz() = wz + F2*gz;
}

void Collisionpair::UpdateParticleVelInfo()
{
    Real sc(sin(chi)), cc(cos(chi));
    Real se(sin(eta)), ce(cos(eta));
    gx = g1 * (ct * cc - st * sc * ce);
    gy = g1 * (stcp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g1 * (stsp * cc + ct * sp * sc * ce + cp * sc * se);
    pt.vx() = wx + F2 * gx;
    pt.vy() = wy + F2 * gy;
    pt.vz() = wz + F2 * gz;
}

void Collisionpair::EjectIonReaction(Particle& particle)
{
    Real vx, vy, vz;
    VelBoltzDistr(vth, vx, vy, vz);
    particle.vx() = vx; 
    particle.vy() = vy; 
    particle.vz() = vz;
}
