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
    Particles::size_type ipart, npart, ncoll;
    
    const int spec_id = (reaction_arr[icsp].first)[0];
    Reaction* & reaction = reaction_arr[icsp].second;
    Particles* & pts = species_arr[spec_id]->particles;
    const Real pm = species_arr[spec_id]->mass;
    const Real m = (pm * mass)/(pm + mass);

    const std::string& name = species_arr[spec_id]->name;
    std::ofstream of(name+".dat", std::ofstream::app);
    std::ofstream coll("coll.dat", std::ofstream::app);
    Real nu_max(0);
    npart = pts->size();

    VecRealArr nu_arr(npart);
    std::vector<VrArr> vr_arr(npart);
    std::vector<Real> vel_r(npart);
    std::vector<int> index_list;
    
    for (ipart=0; ipart < npart; ++ipart) {
        Real nevrt, nutot(0);
        Real vxb, vyb, vzb;
        Real energy;
        Particle& pt = (*pts)[ipart];

        VelBoltzDistr(vth, vxb, vyb, vzb);
        vr_arr[ipart] = {pt.vx()-vxb, pt.vy()-vyb, pt.vz()-vzb};
        vel_r[ipart] = velocity(vr_arr[ipart][0], vr_arr[ipart][1], vr_arr[ipart][2]);
        energy = 0.5 * vel_r[ipart]*vel_r[ipart] * m;
        nu_arr[ipart] = reaction->en_cs(energy);

        nevrt = vel_r[ipart] * ndens;
        for (auto& nui : nu_arr[ipart]) {
            nui *= nevrt;
            nutot += nui;
        }
        if (nutot > nu_max) nu_max = nutot;
    }
    ncoll = static_cast<Particles::size_type>(npart*Pcoll(nu_max,dt)+0.5);
    random_index(npart, ncoll, index_list);
    sort(index_list.begin(), index_list.end());

    CollProd products;
    int ntype = reaction->isize();
    for(const int& ipart: index_list) {
        Particle& ptc = (*pts)[ipart];
        const std::vector<Real>& nu = nu_arr[ipart];
        
        Real rnd = ranf(), nuj = 0.;
        int itype = 0;
        while(itype != ntype) {
            nuj += nu[itype];
            if(rnd < (nuj/nu_max)) {
                Collisionpair collision = Collisionpair(ptc, vr_arr[ipart], vel_r[ipart], pm, mass, vth);
                ParticleCollision(itype, mass, reaction, collision, products);
                // break;
                if (itype == 0) ++ela;
                else if (itype == 1) ++exc;
                else ++ion;
                break;
            }
            ++itype;
        }
    }

    if(!products.empty()){

        std::vector<int>& prodid =  reaction->productid_arr;
        int nprods = static_cast<int>(products.size());
        int nid = static_cast<int>(prodid.size());
        for(int iprod = 0; iprod < nprods; ++iprod) {
            for(int i = 0; i < nid; ++i){
                int ispec = prodid[i];
                species_arr[ispec]->particles->append(products[iprod][i]);
            }
        }
    }
    coll << " nparts: " << npart  << " nu_max: " << nu_max 
         << " ncolls: " << ncoll << " -> "
         << ela << " " << exc << " " << ion << std::endl;

    species_arr[spec_id]->get_particles_energy();
    of << species_arr[spec_id]->toten << std::endl;
    coll.close();
    of.close();
}

void Tile::ParticleColumnCollision(Real dt, int icsp)
{ 
    espic_error("Column collision has not prepared");
}

// void Tile::NullCollisionMethod(Particles& pts, Real dt, Real nu_max, 
//                          Reaction*& react, Real imass, CollProd& products)
// {
//     Particles::size_type ncoll, ipart, npart;
//     pts.particles_shuffle();
//     std::ofstream of("part.out", std::ofstream::app);
    
//     npart = pts.size();
//     ncoll = static_cast<Particles::size_type>(npart*Pcoll(nu_max,dt));
//     std::cout << nu_max << " " << ncoll << " ";
//     int ntype(react->isize());
//     for(ipart = 0; ipart < ncoll; ++ipart) {
//         Particle& ptc = pts[ipart];
//         const std::vector<Real>& nu = ptc.nu();

//         // of << "v = (" << ptc.vx() << " "<< ptc.vy() << " " << ptc.vz() << " ), " 
//         //    << "vr = (" << ptc.vxr() << " "<< ptc.vyr() << " " << ptc.vzr() << " ), "
//         //    << "nu = (";
//         // for (const Real& iu : nu)
//         //     of << iu << " ";
//         // of <<")." << std::endl;
        
//         Real rnd(ranf()), nuj(0.);
//         int itype = 0;
//         while(itype != ntype) {
//             nuj += nu[itype];
//             if(rnd < nuj / nu_max) {
//                 Collisionpair collision = Collisionpair(ptc, imass, mass, vth);
//                 ParticleCollision(itype, mass, react, collision, products);
//                 if(itype==0) ++ela;
//                 else if(itype==1) ++exc;
//                 else ++ion;
//                 break;
//             }
//             ++itype;
//         }
//     }
// }

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
                std::cout << "none ";
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
        std::cout << "]" << std::endl;
    }
    
}
