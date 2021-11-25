#include <iostream>
#include <algorithm>
#include "../espic_info.h"
#include "../espic_math.h"
#include "../mesh.h"
#include "../species.h"
#include "beam.h"

using namespace ESPIC;

/* ---------------- Begin Public Methods ---------------- */

/* ------------------------------------------------------- */

Beam::Beam(int dimension, const BeamDef* const & beamdef, const SpeciesDef* const & specdef) 
    : Inject(beamdef->specid,
             beamdef->ndens,
             beamdef->temp,
             beamdef->vel),
      width_x(beamdef->width_x),
      width_y(beamdef->width_y),
      center{beamdef->center[0], beamdef->center[1], beamdef->center[2]},
      cache_size(0),
      cache_cursor(0)
{
  // vector normal to beam tmat[0]
  Real norm = 0.;
  for (int a = 0; a < 3; a++) norm += beamdef->fnorm[a]*beamdef->fnorm[a];
  norm = sqrt(norm);
  for (int a = 0; a < 3; a++) tmat[0][a] = beamdef->fnorm[a]/norm;

  // thermal velocity
  vth = sqrt(2.*temp/specdef->mass);
  weight = specdef->weight;

  switch (dimension) {
    case 2:
      init_2d();
      ptr_gen_particles = &Beam::gen_particles_2d;
      break;
    case 3:
      init_3d();
      ptr_gen_particles = &Beam::gen_particles_3d;
      break;
    case 5:
      init_axi();
      ptr_gen_particles = &Beam::gen_particles_axi;
      break;
    default:
      espic_error("Simulation must be performed in 2d, 3d or axisymmetric");
  }

// std::cout << "beam: " <<  "species = " << specdef->name
// << ", n = " << ndens
// << ", T = " << temp
// << ", v = [" << vel[0]
// << ", " << vel[1]
// << ", " << vel[2] << "]"
// << ", width = " << width
// << ", center = [" << center[0]
// << ", " << center[1]
// << ", " << center[2] <<"]"
// << ", direction =[" << tmat[0][0]
// << ", " << tmat[0][1]
// << ", " << tmat[0][2] << "]"
// << ", area = " << area
// << ", vth = " << vth
// << ", sn = " << sn
// << ", pn = " << pn
// << ", nflowrate = " << nflowrate << "\n";
 
}

/* ------------------------------------------------------- */

void Beam::gen_particles(Real dt, std::vector<Particle>& particles)
{
  (this->*ptr_gen_particles)(dt, particles);
}

/* ------------------------------------------------------- */

void Beam::print() const
{
  std::cout << "beam: " <<  "species id = " << specid
            << ", n = " << ndens
            << ", T = " << temp
            << ", v = [" << vel[0]
            << ", " << vel[1]
            << ", " << vel[2] << "]"
            << ", width_x = " << width_x
            << ", width_y = " << width_y
            << ", center = [" << center[0]
            << ", " << center[1]
            << ", " << center[2] <<"]"
            << ", direction =[" << tmat[0][0]
            << ", " << tmat[0][1]
            << ", " << tmat[0][2] << "]"
            << ", t1 =[" << tmat[1][0]
            << ", " << tmat[1][1]
            << ", " << tmat[1][2] << "]"
            << ", t2 =[" << tmat[2][0]
            << ", " << tmat[2][1]
            << ", " << tmat[2][2] << "]"
            << ", area = " << area
            << ", vth = " << vth
            << ", sn = " << sn
            << ", pn = " << pn
            << ", nflowrate = " << nflowrate << "\n";
}

// void Beam::inject_2d(Real dt, Particle* particle)
// {
//   Smallint np_gen, ip;
//   Real     np_inj;
//   Real     xnew[3] = { 0., 0., 0.};   // xnew[2] always 0 for 2d
//   Index    icell[3];
// 
//   Particle::OnePart *p = new Particle::OnePart();
//   p->x[0] = p->x[1] = p->x[2] = 0.; // p->x[2] is always 0. for 2d
// 
//   for (auto beam = _beams.begin(); beam != _beams.end(); ++beam) {
//     Real dt_weight = dt/particle->get_one_spec(beam->specid).weight;
//     Vector3 dx = { beam->xhi[0] - beam->xlo[0],
//                    beam->xhi[1] - beam->xlo[1],
//                    beam->xhi[2] - beam->xlo[2] };
//     
//     np_inj = dt_weight*beam->nflowrate;
//     np_gen = (Smallint)(np_inj + beam->nres + ranf());
//     beam->nres += (np_inj - np_gen);
// 
//     p->specid = beam->specid;
//     for (ip = 0; ip < np_gen; ++ip) {
// 
//       // get velocity of a particle
//       getvel2d(beam->vth,
//                beam->sn,
//                beam->pn,
//                beam->direction,
//                beam->nvec,
//                p->v);
// 
//       // get position of a particle
//       for (int c = 0; c < 3; ++c)
//         xnew[c] = beam->xlo[c] + dx[c]*ranf();
//       xnew[beam->direction] += p->v[beam->direction]*dt*ranf();
// 
//       p->cellid = _mesh.find_cell(xnew, icell);
// 
//       p->x[0] = xnew[0];
//       p->x[1] = xnew[1];
// 
//       particle->add_one(p);
//     }   // end for (ip = 0; ip < np_gen; ++ip)
//   }   // end for (auto beam = _beams.begin(); ... )
// 
//   delete p;
// 
//   return;
// }

/* ----------------- End Public Methods ----------------- */

/* ---------------- Begin Private Methods ---------------- */

/* ------------------------------------------------------- */

void Beam::init_2d()
{
  init_plane_2d_axi();

  area = width_x*width_y;  // 2d area (unit length in z)

  precomputed(area);

  return;
}

/* ------------------------------------------------------- */

void Beam::init_axi()
{
  Real pi = ESPIC::PI;

  init_plane_2d_axi();
  
  area = 2.*pi*width_y*center[1];

  precomputed(area);

  return;
}

/* ------------------------------------------------------- */

void Beam::init_3d()
{
}

/* ------------------------------------------------------- */

void Beam::init_plane_2d_axi()
{
  tmat[0][2] = 0.;    // 2d/axi-sym case, only on x-y plane

  // first vector on beam exit surface
  // and perpendicular to the normal vector
  tmat[1][0] = -tmat[0][1];
  tmat[1][1] = +tmat[0][0];
  tmat[1][2] = 0.;

  // vector normal to the two other vectors
  tmat[2][0] = 0.;
  tmat[2][1] = 0.;
  tmat[2][2] = tmat[0][0]*tmat[1][1] - tmat[0][1]*tmat[1][0];
}

/* ------------------------------------------------------- */

void Beam::init_particle_cache_arr(Real dt)
{
  Real dt_weight = dt/weight;
  Real np_inj = dt_weight*nflowrate;

  // cache ~ 1M particles in particle_cache_arr
  cache_size = static_cast<int> (std::min(2*1024.*1024./np_inj, 1000.));
  // round to multiple of 10
  cache_size = static_cast<int> (cache_size/10 + 0.5)*10;
  cache_count.resize(cache_size+1, 0);

  if (np_inj > 1e-3) {  // do not inject particles if np_inj too small 
    for (cache_cursor = 0; cache_cursor < cache_size; cache_cursor++) {
      int np_gen = static_cast<int> (np_inj + nres + ranf());
      cache_count[cache_cursor+1] = cache_count[cache_cursor]+np_gen;
      nres += (np_inj - np_gen);
    }
    particle_cache_arr.resize(cache_count[cache_size]);
  }
}

/* ------------------------------------------------------- */

void Beam::gen_particles_2d(Real dt, std::vector<Particle>& particles)
{
  if (cache_cursor == cache_size) {
    init_particle_cache_arr(dt);

    Real dtfrac = dt*1.0;
    Particle particle (0, 0, 0, 0, 0, 0);
    Real pos[3] = {0., 0., 0.}, vel[3] = {0., 0., 0.}, vn;
    Real Lx, Ly, Lz = 0.;
    Real x0 = center[0] - 0.5*width_y*tmat[0][1]; // origin of beam plane coord system
    Real y0 = center[1] - 0.5*width_y*tmat[0][0];

// #ifdef OMP
// #pragma omp parallel for private (particle, pos,vel, vn, Lx, Ly, Lz)
// #endif
    for (int ip = 0; ip < cache_count[cache_size]; ip++) {
      // generate velocity for a particle to be injected
      gen_one_vel(vel);
      
      // generate initial position for a particle to be injected
      vn = vel[0]*tmat[0][0] + vel[1]*tmat[0][1] + vel[2]*tmat[0][2];
      Lx = vn*ranf()*dtfrac;
      Ly = width_y*ranf();
      pos[0] = tmat[0][0]*Lx + tmat[1][0]*Ly + tmat[2][0]*Lz + x0;
      pos[1] = tmat[0][1]*Lx + tmat[1][1]*Ly + tmat[2][1]*Lz + y0;

#ifdef DEBUG
      if (Lx < 0.) {
        espic_error("Lx cannot be negative for particles to be injected.");
      }
#endif

      particle.x() = pos[0];
      particle.y() = pos[1];
      particle.z() = 0.;
      particle.vx() = vel[0];
      particle.vy() = vel[1];
      particle.vz() = vel[2];

//     particles[ip] = particle;
      particle_cache_arr[ip] = particle;
    }   // end for (ip = 0; ip < np_gen; ++ip)

// if (specid == 0) {
//   for (int ip = 0; ip < cache_count[cache_size]; ip++) {
//     auto particle = particle_cache_arr[ip];
//     std::cout <<  particle.vx() << " " << particle.vy() << " " << particle.x() << " " << particle.y() << "\n";
//   }
// }
    cache_cursor = 0;
  }

  inject_beam_particles(particles);

}

/* ------------------------------------------------------- */

void Beam::gen_particles_axi(Real dt, std::vector<Particle>& particles)
{
  if (cache_cursor == cache_size) {
    init_particle_cache_arr(dt);

    Real dtfrac = dt*1.0;
    Particle particle (0, 0, 0, 0, 0, 0);
    Real pos[3] = {0., 0., 0.}, vel[3] = {0., 0., 0.}, vn;
    Real Lx, Ly, Lz = 0.;
    Real x0 = center[0] - 0.5*width_y*tmat[0][1]; // origin of beam plane coord system
    Real y0 = center[1] - 0.5*width_y*tmat[0][0];
    Real y0sq = y0*y0, y0doub = 2.*y0;
    bool normal_to_x = fabs(tmat[0][0]) < 1e-13;  // beam plane is normal to x-axis
    Real wy_nx = width_y*tmat[0][0];

    for (int ip = 0; ip < cache_count[cache_size]; ip++) {
      // generate velocity for a particle to be injected
      gen_one_vel(vel);

      if (normal_to_x) {
        Ly = width_y*ranf();
      }
      else {
        Ly = -y0 + sqrt(y0sq + wy_nx*(y0doub + wy_nx)*ranf());
      }
      // generate initial position for a particle to be injected
      vn = vel[0]*tmat[0][0] + vel[1]*tmat[0][1] + vel[2]*tmat[0][2];
      Lx = vn*ranf()*dtfrac;
      pos[0] = tmat[0][0]*Lx + tmat[1][0]*Ly + tmat[2][0]*Lz + x0;
      pos[1] = tmat[0][1]*Lx + tmat[1][1]*Ly + tmat[2][1]*Lz + y0;

#ifdef DEBUG
      if (Lx < 0.) {
        espic_error("Lx cannot be negative for particles to be injected.");
      }
#endif

      particle.x() = pos[0];
      particle.y() = pos[1];
      particle.z() = 0.;
      particle.vx() = vel[0];
      particle.vy() = vel[1];
      particle.vz() = vel[2];

      particle_cache_arr[ip] = particle;
    }   // end for (ip = 0; ip < np_gen; ++ip)

    cache_cursor = 0;
  }

  inject_beam_particles(particles);

}


/* ------------------------------------------------------- */

void Beam::gen_particles_3d(Real dt, std::vector<Particle>& particles)
{

}

/* ------------------------------------------------------- */

void Beam::inject_beam_particles(std::vector<Particle>& particles)
{
  auto pos = particles.cbegin();
  auto beg = particle_cache_arr.cbegin()+cache_count[cache_cursor];
  auto end = particle_cache_arr.cbegin()+cache_count[cache_cursor+1];
  particles.insert(pos, beg, end);
  cache_cursor++;
}

/* ----------------- End Private Methods ----------------- */
