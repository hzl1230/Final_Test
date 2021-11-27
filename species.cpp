#include <cstdio>
#include <cstring>

#include "species.h"
#include "espic_info.h"

Species::Species(const SpeciesDef* const & specdef)
  : name (specdef->name),
    mass (specdef->mass),
    charge (specdef->charge),
    weight (specdef->weight),
    particles (new Particles())
{
}

Species::Species(const Species& spec)
  : name(spec.name),
    mass(spec.mass),
    charge(spec.charge),
    weight(spec.weight),
    particles(new Particles(*(spec.particles)))
{
}

Species::~Species()
{
  delete particles;
}

void Species::gen_particles(Bigint n)
{
  std::vector<Real> x(n), vx(n);
  std::vector<Real> y(n), vy(n);
  std::vector<Real> z(n), vz(n);

}

void Species::reserve_num_particles(Bigint n)
{ 
  particles->reserve(n);
}

void Species::write_restart(FILE* fp)
{
  char buf[1024];
  uint32_t size = 0;
  Particles::size_type nparts = num_particles();
  int name_len = static_cast<int> (name.size());
  
  std::memcpy(buf+size, &name_len, sizeof(int));
  size += sizeof(int);
  
  std::memcpy(buf+size, name.c_str(), name_len*sizeof(char));
  size += name_len*sizeof(char);
  
  std::memcpy(buf+size, &mass, sizeof(Real));
  size += sizeof(Real);
 
  std::memcpy(buf+size, &charge, sizeof(Real));
  size += sizeof(Real);
 
  std::memcpy(buf+size, &weight, sizeof(Real));
  size += sizeof(Real);
 
  std::memcpy(buf+size, &nparts, sizeof(Particles::size_type));
  size += sizeof(Particles::size_type);

  fwrite(buf, sizeof(char), size, fp);

  // write particle info
  particles->write_restart(fp);
}
