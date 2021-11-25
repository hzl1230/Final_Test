#include <iostream>

#include "../espic_info.h"
#include "../espic_math.h"
#include "inject.h"
#include "beam.h"

using namespace ESPIC;

// cutoff value for particle injection
const Real vth_max = 3.0;

/* ---------------- class InjectDef ----------------*/ 

/* ------------------------------------------------------- */

InjectDef::InjectDef()
  : beamdef_arr(0)
{

}

/* ------------------------------------------------------- */

InjectDef::~InjectDef()
{
  for (int i = 0; i < num_beams(); i++) delete beamdef_arr[i];
  beamdef_arr.clear();
  beamdef_arr.shrink_to_fit();
}

void InjectDef::append(BeamDef* b)
{
  beamdef_arr.push_back(b);
}

/* ======================================================== */
/* ======================================================== */


/* ---------------- class Inject ----------------*/ 

/* ---------------- Begin Public Methods ---------------- */


/* ---------------- End Public Methods ---------------- */


/* ---------------- Begin Protected Methods ---------------- */

/* ------------------------------------------------------- */

void Inject::precomputed(Real area)
{
  Real h;

  // projection of vd to beam normal
  sn = tmat[0][0]*vel[0] + tmat[0][1]*vel[1] + tmat[0][2]*vel[2];
  sn /= vth;

  h  = sqrt(sn*sn + 2.0);
  pn = 2.0/(sn + h)*exp(0.5 + 0.5*sn*(sn - h));
  nflowrate = area*calc_nflux();

  return;
}

/* ------------------------------------------------------- */

Real Inject::calc_nflux()
{
  Real nflux;

  if (sn >= 0) 
    nflux = 0.5*vth*ndens*(exp(-sn*sn)/sqrt(PI) + sn*(1.0 + erf(sn)));
  else
    nflux = 0.5*vth*ndens*(exp(-sn*sn)/sqrt(PI) + sn*erfc(-sn));

  return nflux;
}

/* ------------------------------------------------------- */

void Inject::gen_one_vel(Real v[3])
{
  Real snew, f, vrf, trf, Lu, Lv, Lw;
  Real vth_range = 2.*vth_max;

  // gen thermal vel (ratio) normal to surface by acceptance-rejection
  while(1) {
    snew = -vth_max + vth_range*ranf();
    if ((snew + sn) <= 0.0) continue;
    f = pn*(snew + sn)*exp(-snew*snew);
    if (f > ranf()) break;
  }

  // tangential thermal vels by Box-Muller method
//   vrf = vth*sqrt(-log(ranf() + SMALLREAL));
  vrf = vth*ranf.normal_dist_factor();
  trf = PI2*ranf();

  // thermal velocity in the surface coordinate system
  Lu = snew*vth;        // normal component
  Lv = vrf*cos(trf);    // tangential components
  Lw = vrf*sin(trf);
    
  // convert v in the surface coordinate system back to the Cartesian system
  v[0] = tmat[0][0]*Lu + tmat[1][0]*Lv + tmat[2][0]*Lw + vel[0];
  v[1] = tmat[0][1]*Lu + tmat[1][1]*Lv + tmat[2][1]*Lw + vel[1];
  v[2] = tmat[0][2]*Lu + tmat[1][2]*Lv + tmat[2][2]*Lw + vel[2];
}

/* ----------------- End Protected Methods ----------------- */