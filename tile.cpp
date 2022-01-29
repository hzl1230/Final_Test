#include "tile.h"
#include "ambient.h"
#include <fstream>

int ela, exc, ion;
using std::cout;
using std::endl;

/*----------------------- begin public method ------------------------*/

/* Constructor */

Tile::Tile(
    Mesh* msh,
    const ParamParticle* param_particle,
    const CrossSection* cross_section)
    : mesh(msh),
      mass(cross_section->background->mass),
      ndens(cross_section->background->ndens),
      vth(cross_section->background->vth),
      ptr_particle_collision(nullptr)
{
    Bigint np = 10000;
    const vector<SpeciesDef*>& specdef_arr = param_particle->specdef_arr;
    const vector<AmbientDef*>& ambdef_arr = param_particle->ambientdef_arr;

    int nspecies = static_cast<int>(specdef_arr.size());
    species_arr.resize(nspecies);
    for (int ispec = 0; ispec < nspecies; ++ispec) {
        species_arr[ispec] = new Species(specdef_arr[ispec]);
        species_arr[ispec]->reserve_num_particles(np);
    }

    InitAmbient(mesh->dimension(), ambdef_arr, specdef_arr);
    InitCollision(param_particle, cross_section);

    if(!ambient_arr.empty()) {
        int num_ambient = static_cast<int>(ambient_arr.size());
        for (int iamb = 0; iamb < num_ambient; ++iamb) {
            Ambient* const& ambient = ambient_arr[iamb];
            Species* & species = species_arr[ambient->species_id()];

            ambient->gen_ambient_0d(np, species->particles);
        }
    }
}

Tile::~Tile()
{
    if(!species_arr.empty()) {
        for (size_t ispec = 0; ispec < species_arr.size(); ++ispec)
            delete species_arr[ispec];
        species_arr.clear();
        species_arr.shrink_to_fit();
    }
    if(!ambient_arr.empty()) { 
        for (size_t iamb = 0; iamb < ambient_arr.size(); ++iamb)
            delete ambient_arr[iamb];
        ambient_arr.clear();
        ambient_arr.shrink_to_fit();
    }
}

void Tile::ParticleCollisioninTiles(Real dt)
{
    size_t num_collspec = reaction_arr.size();
    for (size_t icsp = 0; icsp < num_collspec; ++icsp) {
        std::vector<int>& specid_arr = reaction_arr[icsp].first;
        if (specid_arr.size() < 2) {
            ela = 0; exc = 0; ion = 0;
            ptr_particle_collision = &Tile::ParticleBackgroundCollision;
        }
        else {
            ptr_particle_collision = &Tile::ParticleColumnCollision;
            reaction_arr[icsp].second->is_background_collision = false;
        }
        (this->*ptr_particle_collision)(dt, icsp);
    }
}

void Tile::ParticleBackgroundCollision(Real dt, int icsp)
{
    Particles::size_type npart, ncoll;
    const int spec_id = (reaction_arr[icsp].first)[0];
    Reaction* & reaction = reaction_arr[icsp].second;
    const Real nu_max = reaction->max_coll_freq();
    bool e_incident = reaction->is_e_incident();

    const std::string& name = species_arr[spec_id]->name;
    std::ofstream of(name+".dat", std::ofstream::app);

    Particles* & particles = species_arr[spec_id]->particles;
    std::vector<int> pdep, index_list;

    npart = particles->size();
    ncoll = static_cast<Particles::size_type>(npart*Pcoll(nu_max,dt)+0.5);
    random_index(npart, ncoll, index_list);
    // sort(index_list.begin(), index_list.end());
    for (const Particles::size_type ipart: index_list) {
        reaction->mcc_background_collision(ipart, (*particles)[ipart], species_arr, pdep, e_incident);
    } 
    // RemoveHoles(pdep.size()/2, pdep, particles);
    species_arr[spec_id]->get_particles_energy();
    of << species_arr[spec_id]->toten << std::endl;
    of.close();
}

void Tile::ParticleColumnCollision(Real dt, int icsp)
{ 
    espic_error("Column collision has not prepared");
}



/*----------------------- begin private method ------------------------*/

void Tile::InitAmbient(
    int dimension,
    const std::vector<AmbientDef*>& ambdef_arr,
    const std::vector<SpeciesDef*>& specdef_arr)
{
    int num_ambient = static_cast<int>(ambdef_arr.size());

    for (int i = 0; i < num_ambient; i++) {
        const AmbientDef* const& ambdef = ambdef_arr[i];
        const SpeciesDef* const& specdef = specdef_arr[ambdef->specid];
        ambient_arr.push_back(new Ambient(dimension, ambdef, specdef));
    }
}

void Tile::InitCollision(
        const ParamParticle* pp,
        const CrossSection* cs)
{
    for (int icsp = 0; icsp < cs->num_pairs(); ++icsp) {
        Reaction* reaction = cs->reaction_arr[icsp];
        const ReactPair& spair = reaction->pair();
        // const StringList& prod_list = cs->product_arr[icsp];
        std::vector<int> spec_id, prod_id;
        int specid1 = -1, specid2 = -1;
        Real m1, m2, nref = ndens;
        try {
            specid1 = pp->map_spec_name_indx.at(spair.first);
            spec_id.push_back(specid1);
            m1 = pp->specdef_arr[specid1]->mass;
            if (spair.second != "default"){
                specid2 = pp->map_spec_name_indx.at(spair.second);
                spec_id.push_back(specid2);
                m2 = pp->specdef_arr[specid2]->mass; 
                nref = 1.0;
            } else { m2 = mass; }
            reaction->mr() = (reaction->is_e_incident()) ? m1 : m1*m2 / (m1+m2);
        }
        catch (const std::out_of_range& oor) {
            std::ostringstream oss;
            oss << "Unknown species \"" << spair.first << " " << spair.second
                << "\" given to \"cross_section\" command in [csection.in]";
            espic_error(oss.str());
        }
         // pbc directly get max nu;
         // ppc get max vsigma, n multiply in tile
        reaction->init_collision(m1, m2, nref);
        
        reaction_arr.push_back(std::make_pair(spec_id, reaction));
        std::cout << "Reaction " << icsp  <<", relative mass: " << reaction->mr()
                  << ", Max Coll Freq: " << reaction->max_coll_freq()
                  << " product(name,specid): [";
        
        int rnum = reaction->isize();
        reaction->prodid_arr.resize(rnum);
        for (int irct = 0; irct < rnum; ++irct) {
            const StringList& tempprod = reaction->get_prod(irct);
            std::cout << " ";
            if (tempprod.empty()) {
                std::cout << "none";
                continue;
            }
            for (const std::string& pro : tempprod) {
                int spid = -1;
                try {
                    spid = pp->map_spec_name_indx.at(pro);
                }
                catch (const std::out_of_range& oor) {
                    std::ostringstream oss;
                    oss << "Unknown species \"" << pro << " " << spair.second
                        << "\" given to \"cross_section\" command in [csection.in]";
                    espic_error(oss.str());
                }
                reaction->prodid_arr[irct].emplace_back(spid);
                std::cout << pro << "(" << spid << ")" ;
            }
            
        }  
        std::cout << " ]" << std::endl;
    }
    
}
