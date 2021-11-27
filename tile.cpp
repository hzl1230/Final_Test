#include "tile.h"
#include "ambient.h"

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
        if (specid_arr.size() < 2) 
            ptr_particle_collision = &Tile::ParticleBackgroundCollision;
        else {
            ptr_particle_collision = &Tile::ParticleColumnCollision;
            reaction_arr[icsp].second->is_background_collision = false;
        }
        (this->*ptr_particle_collision)(dt, icsp);
    }
}

void Tile::ParticleBackgroundCollision(Real dt, int icsp)
{
    Particles::size_type ipart, npart;

#ifdef OMP
#pragma omp parallel for private(ipart, npart) \
  schedule(dynamic)
#endif

    const int spec_id = (reaction_arr[icsp].first)[0];
    Reaction* & reaction = reaction_arr[icsp].second;
    Particles* & pts = species_arr[spec_id]->particles;
    const Real pm = species_arr[spec_id]->mass;
    // const std::string& name = species_arr[spec_id]->name;
    Real nu_max(0);
    
    npart = pts->size();
    for (ipart=0; ipart < npart; ++ipart) {
        Real vxb, vyb, vzb;
        Real vel, energy;
        Real nevrt, nutot(0);
        Particle& pt = (*pts)[ipart];
        std::vector<Real> info;
        std::vector<Real>& nu = pt.nu();
        if(!info.empty()) info.clear();
        if(!nu.empty()) nu.clear();

        VelBoltzDistr(vth, vxb, vyb, vzb);
        RelativeVelocity(pt, vxb, vyb, vzb);
    
        energy = get_energy(pt.vxr(), pt.vyr(), pt.vzr());
        info = reaction->en_cs(energy);
        vel = sqrt(2.0 * energy);
        nevrt = vel*ndens;

        transform(info.begin(),info.end(), back_inserter(nu), 
                    [=](auto& x) { return x*nevrt; });
        for(auto& nui : nu) nutot += nui;
        if (nutot > nu_max) nu_max = nutot;
    }
        CollProd products;
        NullCollisionMethod(*pts, dt, nu_max, reaction, pm, products);
        if(!products.empty()){
            for(auto iprod = 0; iprod<products.size(); ++iprod) {
                pts->append(products[iprod][0]);
                species_arr[spec_id+1]->particles->append(products[iprod][1]);
            }
        }
}

void Tile::ParticleColumnCollision(Real dt, int icsp)
{ 
    espic_error("Column collision has not prepared");
}

void Tile::NullCollisionMethod(Particles& pts, Real dt, Real nu_max, 
                         Reaction*& react, Real imass, CollProd& products)
{
    Particles::size_type ncoll, ipart, npart;
    npart = pts.size();
    ncoll = static_cast<Particles::size_type>(npart*Pcoll(nu_max,dt));
    int ntype(react->isize());
    for(ipart = 0; ipart < ncoll; ++ipart) {
        Particle& ptc = pts[ipart];
        const std::vector<Real>& nu = ptc.nu();
        Real rnd(ranf()), nuj(0.);
        int itype(0);
        while(itype != ntype) {
            nuj += nu[itype];
            if(rnd < nuj / nu_max) {
                Collisionpair collision = Collisionpair(ptc, imass, mass, vth);
                ParticleCollision(itype, mass, react, collision, products);
                ++itype;  break;}
        }
    }
}

void Tile::ParticleCollision(
    const int type_id, 
    Real mass,
    Reaction*& reaction,
    Collisionpair& cop,
    CollProd& prod_arr)
{
    Real threshold;
    std::string type = (reaction->get_types())[type_id];
    bool is_bc = reaction->is_background_collision;
    if (type_id == 0)  threshold = 0.0;
    else  threshold = (reaction->th())[type_id-1];

    if ("ela" == type)
        cop.ParticleElasticCollision();
    else if ("exc" == type)
        cop.ParticleExcitatinCollision(threshold);
    else if ("ion" == type) {
        if (is_bc)
          cop.ParticleIonizationCollision(threshold);
        prod_arr.emplace_back(cop.ion_products()); 
    }
    else if ("iso" == type) 
        cop.ParticleIsotropicCollision();
    else if ("back" == type)
        cop.ParticleBackwardCollision();
    else
        espic_error("Unknown Collision Type");
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
        const StringList& prod_list = cs->product_arr[icsp];
        std::vector<int> spec_id, prod_id;
        int specid1 = -1, specid2 = -1;
        try {
            specid1 = pp->map_spec_name_indx.at(spair.first);
            spec_id.emplace_back(specid1);
            if (spair.second != "default"){
                specid2 = pp->map_spec_name_indx.at(spair.second);
                spec_id.emplace_back(specid2);
            }
        }
        catch (const std::out_of_range& oor) {
            std::ostringstream oss;
            oss << "Unknown species \"" << spair.first << " " << spair.second
                << "\" given to \"cross_section\" command in [csection.in]";
            espic_error(oss.str());
        }
        reaction_arr.emplace_back(std::make_pair(spec_id, reaction));
        std::cout << "Reaction " << icsp << " product(name,specid): [ ";
        
        for (const auto& pro: prod_list) {
            if ("none" == pro) {
                std::cout << "none";
                break;
            }
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
            reaction->productid_arr.emplace_back(spid);
            std::cout << pro << "(" << spid << ")" << " ";
        }  
        std::cout << " ]" << std::endl;
    }
    
}
