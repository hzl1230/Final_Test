#include <iostream>
#include <cmath>

#include "espic_info.h"
#include "espic_math.h"
#include "ambient.h"
#include "species.h"
#include "particles.h"

ESPIC::Random ranf;

using namespace ESPIC;

/* ---------------- Begin Public Methods ---------------- */

/* Constructor */
Ambient::Ambient(int dimension, const AmbientDef* const & ambdef, const SpeciesDef* const & specdef)
  : specid (ambdef->specid),
    ndens (ambdef->ndens),
    vth (sqrt(2.*ambdef->temp/specdef->mass)),
    vel {ambdef->vel[0], ambdef->vel[1], ambdef->vel[2]},
    weight (specdef->weight),
    bound_lo {ambdef->bound_lo[0], ambdef->bound_lo[1], ambdef->bound_lo[2]},
    bound_hi {ambdef->bound_hi[0], ambdef->bound_hi[1], ambdef->bound_hi[2]},
    ptr_gen_ambient(nullptr)
{
  switch (dimension) {
    case 2:
      ptr_gen_ambient = &Ambient::gen_ambient_2d;
      break;
    case 3:
      ptr_gen_ambient = &Ambient::gen_ambient_3d;
      break;
    case 5:
      ptr_gen_ambient = &Ambient::gen_ambient_axi;
      break;
    default:
      espic_error("Simulation must be performed in 2d, 3d or axisymmetric");
  }
}

/* ------------------------------------------------------- */

void Ambient::gen_ambient(int nc[3], Real bound_lo[3], Real bound_hi[3], Real dx[3], Particles* &particles)
{
  (this->*ptr_gen_ambient)(nc, bound_lo, bound_hi, dx, particles);
}

void Ambient::gen_ambient_0d(Bigint np, Particles* &particles)
{
  // Bigint np = 10000;
  Particles::size_type nparts, ipart;
  Real x_dim = (bound_hi[0]-bound_lo[0]), y_dim = (bound_hi[1]-bound_lo[1]);
  std::vector<Real> x(np), vx(np);
  std::vector<Real> y(np), vy(np);
  std::vector<Real> z(np), vz(np);
  Real vrf, trf;

  nparts = static_cast<Particles::size_type>(np);
  for (ipart = 0; ipart < nparts; ++ipart){
    x[ipart] = bound_lo[0] + ranf()*x_dim;
    y[ipart] = bound_lo[1] + ranf()*y_dim;
    z[ipart] = 0.;
  }

  for (ipart = 0; ipart < nparts; ipart++) {
    vrf = vth*ranf.normal_dist_factor();
    trf = PI2*ranf();
    vx[ipart] = vrf*cos(trf) + vel[0];
    vy[ipart] = vrf*sin(trf) + vel[1];

    vrf = vth * ranf.normal_dist_factor();
    trf = PI2*ranf();
    vz[ipart] = vrf*cos(trf) + vel[2];
  }
  particles->append(Particles(x, y, z, vx, vy, vz));
}



/* ---------------- End Public Methods ---------------- */

/* ---------------- Begin Private Methods ---------------- */

/* ------------------------------------------------------- */

void Ambient::gen_ambient_2d(int nc[3], Real bound_lo[3], Real bound_hi[3], Real dx[3], Particles* &particles)
{
  Real fparts;
  Particles::size_type nparts, ipart;
  Real x_dim = (bound_hi[0] - bound_lo[0]), y_dim = (bound_hi[1]-bound_lo[1]);
  Real vrf, trf;

  fparts = ndens*x_dim*y_dim/weight;
  if (fparts <= 0.) return;

  nparts = static_cast<Particles::size_type>(fparts + ranf());
  std::vector<Real> x(nparts), y(nparts), z(nparts), vx(nparts), vy(nparts), vz(nparts);
  for (ipart = 0; ipart < nparts; ipart++) {
    x[ipart] = bound_lo[0] + ranf()*x_dim;
    y[ipart] = bound_lo[1] + ranf()*y_dim;
  }

  for (ipart = 0; ipart < nparts; ipart++) {
//     vrf = vth*sqrt(-log(ranf() + SMALLREAL));
    vrf = vth*ranf.normal_dist_factor();
    trf = PI2*ranf();
    vx[ipart] = vrf*cos(trf) + vel[0];
    vy[ipart] = vrf*sin(trf) + vel[1];
  }

  particles->append(Particles(x, y, z, vx, vy, vz));
}

/* ------------------------------------------------------- */

void Ambient::gen_ambient_3d(int nc[3], Real bound_lo[3], Real bound_hi[3], Real dx[3], Particles* &particles)
{
}

/* ------------------------------------------------------- */

void Ambient::gen_ambient_axi(int nc[3], Real bound_lo[3], Real bound_hi[3], Real dx[3], Particles* &particles)
{
  Real fparts;
  Particles::size_type nparts, ipart;
  Real x_dim = (bound_hi[0] - bound_lo[0]);
  Real r0sq = bound_lo[1]*bound_lo[1];
  Real r1sq = bound_hi[1]*bound_hi[1];
  Real drsq = r1sq - r0sq;
  Real vrf, trf;

  fparts = PI*ndens*x_dim*drsq/weight;
  if (fparts <= 0.) return;

  nparts = static_cast<Particles::size_type>(fparts + ranf());
  std::vector<Real> x(nparts), y(nparts), z(nparts, 0), vx(nparts), vy(nparts), vz(nparts, 0);
  for (ipart = 0; ipart < nparts; ipart++) {
    x[ipart] = bound_lo[0] + ranf()*x_dim;
    y[ipart] = sqrt(r0sq + drsq*ranf());
  }

  for (ipart = 0; ipart < nparts; ipart++) {
//     vrf = vth*sqrt(-log(ranf() + SMALLREAL));
    vrf = vth*ranf.normal_dist_factor();
    trf = PI2*ranf();
    vx[ipart] = vrf*cos(trf) + vel[0];
    vy[ipart] = vrf*sin(trf) + vel[1];

//     vrf = vth*sqrt(-log(ranf() + SMALLREAL));
    vrf = vth*ranf.normal_dist_factor();
    trf = PI2*ranf();
    vz[ipart] = vrf*cos(trf) + vel[2];
  }

  particles->append(Particles(x, y, z, vx, vy, vz));
}

