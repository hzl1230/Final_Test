#ifndef _SPECIES_H
#define _SPECIES_H

#include <vector>
#include <string>

#include "espic_type.h"
#include "particles.h"

class SpeciesDef {
  public:
    // Constructor
    SpeciesDef (const std::string& nm, Real m, Real q, Real w)
      : name(nm),
        mass(m),
        charge(q),
        weight(w) {}

    // Copy constructor
    SpeciesDef (const SpeciesDef& orig)
      : name(orig.name),
        mass(orig.mass),
        charge(orig.charge),
        weight(orig.weight) {}

    std::string name;
    Real mass;
    Real charge;
    Real weight;
};

class Species {
  public:
    /* Constructor */
    explicit Species(const SpeciesDef* const &);

    /* Copy constructor */
    Species(const Species&);

    /* Destructor */
    ~Species();

    /* Public methods */
    void reserve_num_particles(Bigint n);

    Particles::size_type num_particles() const { return particles->size(); }

    void write_restart(FILE *);

    /* Members */
    std::string name;
    Real mass;
    Real charge;
    Real weight;
    class Particles* particles;
};

#endif
