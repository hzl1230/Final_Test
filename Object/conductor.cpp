
#include "conductor.h"

int Object::current_id = 0;

/* ---------------- Begin Public Methods ---------------- */

void Conductor::lostparticles_increment(int nspecies,
                                        Real curr_time)
{
  if (lostparticles.back().is_used()) {
    lostparticles.push_back(LostParticleInfo(nspecies, curr_time));
  }
  else  {
    lostparticles.back().time = curr_time;
  }
}

/* ------------------------------------------------------- */

bool Conductor::scrape_particle(int spec_id,
                                Real curr_time,
                                const Vector3& pos_old,
                                const Vector3& pos_new)
{
  bool is_scraped = intersect(pos_old, pos_new);
  if (is_lost_particles_statistics_on()) {
    if (is_scraped) {
      lostparticles.back().add_one_lost(spec_id);
    }
  }
  return is_scraped;
}

/* ----------------- End Public Methods ----------------- */
