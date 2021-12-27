#ifndef _TILE_H_
#define _TILE_H_
#include "cross_section.h"
#include "species.h"
#include "param_particle.h"
#include "collision.h"
#include "mesh.h"

typedef size_t size_type; 
using std::vector;
using std::string;
using std::pair;

class Tile {
public:
    Tile(class Mesh*,
         const class ParamParticle*,
         const class CrossSection*);
    
    ~Tile();

    void ParticleCollisioninTiles(Real);

    void ParticleBackgroundCollision(Real dt, int icps);

    void ParticleColumnCollision(Real dt, int icps);

    void NullCollisionMethod(Particles&, Real, Real, Reaction*&, Real, CollProd&);

    void ParticleCollision(const int , Real ,
                           Reaction*& ,
                           Collisionpair& ,
                           CollProd&);

    typedef std::array<Real, 3> VrArr;
    typedef std::vector<std::vector<Real>> VecRealArr;

private:

    void InitAmbient(int, 
        const vector<AmbientDef*>&, 
        const vector<SpeciesDef*>&);

    void InitCollision(
         const class ParamParticle*,
         const class CrossSection*);
    

    class Mesh* mesh;
    vector<class Ambient*> ambient_arr;
    vector<pair<vector<int>, class Reaction*>> reaction_arr;
    vector<class Species*> species_arr;
    const Real mass, ndens, vth;
    Real xmin, ymin, zmin, xmax, ymax, zmax;
    Real dx, dy, dz;
    Real dxinv, dyinv, dzinv;

    typedef void (Tile::*ParticleCollisioninTile)(Real, int);
    ParticleCollisioninTile ptr_particle_collision;
};

#endif