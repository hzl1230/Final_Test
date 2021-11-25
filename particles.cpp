#include <iostream> 
#include <cstdio>
#include <cstring>
#include <random>
#include <algorithm>

#include "espic_info.h"
#include "particles.h"

/* ---------------- Begin Public Methods ---------------- */

/* Constructors */

/* ------------------------------------------------------- */

Particles::Particles()
  : nparticles(0),
    data(0),
    scalar(0)
//     pos_x(0),
//     pos_y(0),
//     pos_z(0),
//     vel_x(0),
//     vel_y(0),
//     vel_z(0)
{
}

/* ------------------------------------------------------- */

Particles::Particles(const std::vector<Real>& x,
                     const std::vector<Real>& y,
                     const std::vector<Real>& z,
                     const std::vector<Real>& vx,
                     const std::vector<Real>& vy,
                     const std::vector<Real>& vz)
  : nparticles(x.size())
//     pos_x(x),
//     pos_y(y),
//     pos_z(z),
//     vel_x(vx),
//     vel_y(vy),
//     vel_z(vz)
{
  if (nparticles != y.size())
    espic_error("Failed to construct particles because list length is not matched");
  if (nparticles != z.size())
    espic_error("Failed to construct particles because list length is not matched");
  if (nparticles != vx.size())
    espic_error("Failed to construct particles because list length is not matched");
  if (nparticles != vy.size())
    espic_error("Failed to construct particles because list length is not matched");
  if (nparticles != vz.size())
    espic_error("Failed to construct particles because list length is not matched");

  data.resize(nparticles);
  for (size_type ip = 0; ip < nparticles; ip++) {
    data[ip].x() = x[ip];
    data[ip].y() = y[ip];
    data[ip].z() = z[ip];
    data[ip].vx() = vx[ip];
    data[ip].vy() = vy[ip];
    data[ip].vz() = vz[ip];
  }
}

/* Copy Constructor */

/* ------------------------------------------------------- */

Particles::Particles(const Particles& other)
  : nparticles(other.nparticles),
    data(other.data)
//     pos_x(other.pos_x),
//     pos_y(other.pos_y),
//     pos_z(other.pos_z),
//     vel_x(other.vel_x),
//     vel_y(other.vel_y),
//     vel_z(other.vel_z)
{
}

/* Destructor */

/* ------------------------------------------------------- */

Particles::~Particles()
{
//   std::cout << "In ~Particles(): nparticles = " << nparticles << std::endl;
  nparticles = 0;
}

/* ------------------------------------------------------- */

void Particles::particles_shuffle() 
    { 
        srand((unsigned int)time(NULL));
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(data.begin(), data.end(), g);
    }

/* ------------------------------------------------------- */

void Particles::get_sub_particles(size_type n, Particles& sub)
    {
        particles_shuffle();
        if (sub.size() != 0) {
            sub.data.clear();
            sub.nparticles = 0;
        }   
        std::copy(data.begin(), data.begin() + n, back_inserter(sub.data));
        sub.nparticles += n;
    }

/* ------------------------------------------------------- */

void Particles::reserve(size_type n)
{
//   pos_x.reserve(n);
//   pos_y.reserve(n);
//   pos_z.reserve(n);
//   vel_x.reserve(n);
//   vel_y.reserve(n);
//   vel_z.reserve(n);
  data.reserve(n);
}

/* ------------------------------------------------------- */

void Particles::append(const Particle& p)
{
//   pos_x.push_back(p.x());
//   pos_y.push_back(p.y());
//   pos_z.push_back(p.z());
//   vel_x.push_back(p.vx());
//   vel_y.push_back(p.vy());
//   vel_z.push_back(p.vz());
  data.push_back(p);
  ++nparticles;
}

/* ------------------------------------------------------- */

void Particles::append(const std::vector<Particle>& p_arr)
{
  for (size_type i = 0; i < p_arr.size(); i++) append(p_arr[i]);
  nparticles += p_arr.size();
}

/* ------------------------------------------------------- */

void Particles::append(std::vector<Particle>::const_iterator beg,
                       std::vector<Particle>::const_iterator end)
{
  for (auto it = beg; it != end; ++it) {
//     pos_x.push_back(it->x());
//     pos_y.push_back(it->y());
//     pos_z.push_back(it->z());
//     vel_x.push_back(it->vx());
//     vel_y.push_back(it->vy());
//     vel_z.push_back(it->vz());
    data.push_back(*it);
    nparticles++;
  }
}

/* ------------------------------------------------------- */

void Particles::append(std::vector<Particle>::iterator beg,
                       std::vector<Particle>::iterator end)
{
  for (auto it = beg; it != end; ++it) {
//     pos_x.push_back(it->x());
//     pos_y.push_back(it->y());
//     pos_z.push_back(it->z());
//     vel_x.push_back(it->vx());
//     vel_y.push_back(it->vy());
//     vel_z.push_back(it->vz());
    data.push_back(*it);
    nparticles++;
  }
}

/* ------------------------------------------------------- */

void Particles::append(const Particles& others)
{
//   pos_x.insert(pos_x.end(), others.pos_x.begin(), others.pos_x.end());
//   pos_y.insert(pos_y.end(), others.pos_y.begin(), others.pos_y.end());
//   pos_z.insert(pos_z.end(), others.pos_z.begin(), others.pos_z.end());
//   vel_x.insert(vel_x.end(), others.vel_x.begin(), others.vel_x.end());
//   vel_y.insert(vel_y.end(), others.vel_y.begin(), others.vel_y.end());
//   vel_z.insert(vel_z.end(), others.vel_z.begin(), others.vel_z.end());

  data.insert(data.end(), others.data.begin(), others.data.end());
  nparticles += others.size();
}

/* ------------------------------------------------------- */

void Particles::append(Particles* particles)
{
  append(*particles);
}

void Particles::erase(size_type id)
{
//   pos_x.at(id) = pos_x.back();
//   pos_y.at(id) = pos_y.back();
//   pos_z.at(id) = pos_z.back();
//   vel_x.at(id) = vel_x.back();
//   vel_y.at(id) = vel_y.back();
//   vel_z.at(id) = vel_z.back();
//   pos_x.pop_back();
//   pos_y.pop_back();
//   pos_z.pop_back();
//   vel_x.pop_back();
//   vel_y.pop_back();
//   vel_z.pop_back();
  data[id] = data.back();
  data.pop_back();
  nparticles--;
}

/* ------------------------------------------------------- */

void Particles::erase(size_type id, size_type n)
{
  n = id + n <= nparticles ? n : nparticles - id;

  for (size_type i = std::max(id+n, nparticles-n); i < nparticles; i++, id++) {
//     pos_x.at(id) = pos_x.at(i);
//     pos_y.at(id) = pos_y.at(i);
//     pos_z.at(id) = pos_z.at(i);
//     vel_x.at(id) = vel_x.at(i);
//     vel_y.at(id) = vel_y.at(i);
//     vel_z.at(id) = vel_z.at(i);
    data[id] = data[i];
  }
  pop_back(n);
}

/* ------------------------------------------------------- */

void Particles::pop_back(size_type n)
{
//   pos_x.erase(pos_x.end()-n, pos_x.end());
//   pos_y.erase(pos_y.end()-n, pos_y.end());
//   pos_z.erase(pos_z.end()-n, pos_z.end());
//   vel_x.erase(vel_x.end()-n, vel_x.end());
//   vel_y.erase(vel_y.end()-n, vel_y.end());
//   vel_z.erase(vel_z.end()-n, vel_z.end());
  data.erase(data.end()-n, data.end());
  nparticles -= n;
}

/* ------------------------------------------------------- */

void Particles::write_restart(FILE* fp)
{
//   std::vector<Real>::pointer ptr;
//   size_type off = 0, size_bytes = nparticles*sizeof(Real);
//   constexpr uint16_t nprops = 6;
//   char *buf = new char [nprops*size_bytes];
// 
//   ptr = pos_x.data();
//   std::memcpy(buf+off, ptr, size_bytes);
//   off += size_bytes;
// 
//   ptr = pos_y.data();
//   std::memcpy(buf+off, ptr, size_bytes);
//   off += size_bytes;
// 
//   ptr = pos_z.data();
//   std::memcpy(buf+off, ptr, size_bytes);
//   off += size_bytes;
// 
//   ptr = vel_x.data();
//   std::memcpy(buf+off, ptr, size_bytes);
//   off += size_bytes;
// 
//   ptr = vel_y.data();
//   std::memcpy(buf+off, ptr, size_bytes);
//   off += size_bytes;
// 
//   ptr = vel_z.data();
//   std::memcpy(buf+off, ptr, size_bytes);
//   off += size_bytes;

  std::vector<Particle>::pointer ptr;
  size_type size_bytes = nparticles*sizeof(Particle);
  char *buf = new char [size_bytes];

  ptr = data.data();
  std::memcpy(buf, ptr, size_bytes);
  fwrite(buf, size_bytes, 1, fp);

  delete [] buf;
}

/* ------------------------------------------------------- */

void Particles::read_restart_and_store(size_type np, FILE* fp)
{
//   std::vector<Real>::pointer ptr;
//   size_type off = 0, size_bytes = np*sizeof(Real);
//   constexpr uint16_t nprops = 6;
//   char *buf = new char [nprops*size_bytes];
// 
//   size_t sizeread = fread(buf, sizeof(char), nprops*size_bytes, fp);
//   if (sizeread != nprops*size_bytes) {
//     espic_error("Failed to read restart");
//   }
// 
//   pos_x.resize(np);
//   ptr = pos_x.data();
//   std::memcpy(ptr, buf+off, size_bytes);
//   off += size_bytes;
// 
//   pos_y.resize(np);
//   ptr = pos_y.data();
//   std::memcpy(ptr, buf+off, size_bytes);
//   off += size_bytes;
// 
//   pos_z.resize(np);
//   ptr = pos_z.data();
//   std::memcpy(ptr, buf+off, size_bytes);
//   off += size_bytes;
// 
//   vel_x.resize(np);
//   ptr = vel_x.data();
//   std::memcpy(ptr, buf+off, size_bytes);
//   off += size_bytes;
// 
//   vel_y.resize(np);
//   ptr = vel_y.data();
//   std::memcpy(ptr, buf+off, size_bytes);
//   off += size_bytes;
// 
//   vel_z.resize(np);
//   ptr = vel_z.data();
//   std::memcpy(ptr, buf+off, size_bytes);
//   off += size_bytes;

  std::vector<Particle>::pointer ptr;
  const size_type size_bytes = np*sizeof(Particle);
  char *buf = new char [size_bytes];

  size_t sizeread = fread(buf, sizeof(char), size_bytes, fp);
  if (sizeread != size_bytes) {
    espic_error("Failed to read restart");
  }

  data.resize(np);
  ptr = data.data();
  std::memcpy(ptr, buf, size_bytes);

  nparticles = np;

  delete [] buf;
}
/* ---------------- End Public Methods ---------------- */
