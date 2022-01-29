#include "reaction.h"

using std::cout;
using std::endl;

Reaction::Reaction (std::string file, const ReactPair& spair, const Real vtb, int id) 
: infile(file),
spec_pair(spair),
reaction_id(id),
threshold(0),
vth(vtb), mr_(0), 
nu_max(0), 
is_e_in(false),
collision(nullptr)
{
    FILE* fp = fopen(file.c_str(), "r");
    std::vector<std::string> line;
    if (NULL == fp) {
        std::ostringstream oss;
        oss << "Cannot read file [" << infile << "]";
        espic_error(oss.str());
    }
    while(ParseLine(line, fp)){
        if(!info.empty()) info.clear();
        if (line.empty()) continue; 
        else if("basic"==line.at(0)) {
            if (line.size() > 4) {
                info_size = atoi(line[1].c_str());
                arr_length = atoi(line[2].c_str()); 
                de_ = static_cast<Real>(atof(line[3].c_str()));
                deinv_ = 1/de_;
                n_sub = atoi(line[4].c_str());
                element_resize();
            } else {
                std::ostringstream oss;
                oss << "Command basic in file [" << infile << "] need more parameters" ;
                espic_error(oss.str());
            }
        }
        else if("reaction"==line.at(0)) {
            int id = atoi(line[1].c_str());
            if (id > info_size) espic_error("Too Many Reaction");
            types[id-1] = line[2];
            threshold[id-1] = static_cast<Real>(atof(line[3].c_str()));
            line.erase(line.begin(),line.begin()+4);
            if (!line.empty()) {
                int rest = static_cast<int>(line.size());
                for (int i = 0; i < rest; ++i) 
                    product_arr[id-1].emplace_back(line[i]);
            }
        }
        else{
            energy.emplace_back((Real)atof(line[0].c_str()));
            transform(line.begin()+1, line.begin()+info_size+1, back_inserter(info), toReal);
            info_arr.emplace_back(std::move(info));
        }
    }
    cout << "Initial Particle Reaction " << reaction_id << ": " << endl;
    cout << "Reactant: [" << spec_pair.first << "," << spec_pair.second << "],\n"
            << " Threshold: [ " ;
    for (const auto& th: threshold) 
        cout << th << " ";
    cout << "]. " ;
    cout << "de: " << de_ << ", " << "Cross Section Number: " 
            << info_size << ".\n Reaction Type: [";
    for (const auto& type: types)
        cout << " " << type ;
    cout << " ]" << endl;
    cout << "Process collision in every " << n_sub << " steps." << endl;

    if ("e" == spair.first) is_e_in = true;
    is_background_collision = true;
}    

/* ------------------------------------------------------------------------- */

Reaction::~Reaction()
{ 
    if (nullptr != collision) {
        delete collision;
        collision = nullptr;
    }
}

/* ------------------------------------------------------------------------- */

void Reaction::init_collision(const Real m1, const Real m2, const Real n)
{
    collision = new Collisionpair(m1, m2, vth);
    find_max_coll_freq(n);
}

/* -------------------------------------------------------------------------*/

void Reaction::mcc_background_collision(const int ipart,
                                        Particle& particle, 
                                        std::vector<Species*>& species_arr,
                                        std::vector<int>& pdep,
                                        bool e_incident)
{
    std::vector<Real> nu(info_size);
    Real* v = particle.vel();
    Real vel_R;
    if (!is_e_in) {
        Real vxb, vyb, vzb;
        VelBoltzDistr(vth, vxb, vyb, vzb);
        v[0] -= vxb; v[1] -= vyb; v[2] -= vzb;
    }
    vel_R = velocity(v[0], v[1], v[2]);
    const Real energy = 0.5 * vel_R*vel_R * mr_;
    CrossSectionInterplot(energy, nu);
    Real rnd = ranf(), nuj = 0.;
    for (int itype = 0; itype < info_size; itype++) {
        nuj += nu[itype];
        if(rnd < (nuj/nu_max)) {
            std::string type = types[itype];
            const Real th = threshold[itype];
            collision->gen_collision_info(particle, v, vel_R, energy, e_incident);
            if ("ela" == type)
                collision->ElasticScatter(particle);
            else if ("exc" == type || "vib" == type || "rot" == type)
                collision->ExcitatinCollision(th, particle);
            else if ("ion" == type) 
                collision->IonizationCollision(th, particle, prodid_arr[itype], species_arr);
            else if ("disa" == type) { 
                // Generate new particles and remove incident particle 
                collision->IonizationCollision(th, particle, prodid_arr[itype], species_arr);
                pdep.push_back(ipart);
                pdep.push_back(9);
            }
            else if ("det" == type) {
                collision->DetachmentCollision(th, particle, prodid_arr[itype], species_arr);
                pdep.push_back(ipart);
                pdep.push_back(9);
            }
            else if ("iso" == type) 
                collision->IsotropicCollision(particle);
            else if ("back" == type)
                collision->BackwardCollision(particle);
            else if ("cex" == type)
                collision->ChargeExchange(particle);
            else
                espic_error("Unknown Collision Type in Process P-B Collision");
    
            break;
        }
    }

}
void Reaction::mcc_particle_collision(const int ipart, 
                                      const int tpart,
                                      Particle& ipt, 
                                      Particle& tpt,
                                      std::vector<Species*>& species_arr, 
                                      std::vector<int>& ipdep,
                                      std::vector<int>& tpdep)
{
    std::vector<Real> nu(info_size);
    Real v[3] = {ipt.vx()-tpt.vx(), ipt.vy()-tpt.vy(), ipt.vz()-tpt.vz()};
    Real vel_R = velocity(v[0], v[1], v[2]);
    const Real energy = 0.5 * vel_R*vel_R * mr_;
    CrossSectionInterplot(energy, nu);
    Real rnd = ranf(), nuj = 0.0;
    for (int itype = 0; itype < info_size; itype++) {
        nuj += nu[itype];
        if(rnd < (nuj/nu_max)) {
            std::string type = types[itype];
            const Real th = threshold[itype];
            collision->gen_collision_info(ipt, v, vel_R, energy, 0);
            if ("ela" == type)
                collision->CoulombScatter(ipt, tpt);
            else if ("drc" == type || "mn" == type) {
                ipdep.push_back(ipart); ipdep.push_back(9);
                tpdep.push_back(tpart); tpdep.push_back(9);
            }
            else if ("eid" == type) {
                collision->ElectronImpactDetechment(th, ipt, prodid_arr[itype], species_arr);
                tpdep.push_back(tpart);
                tpdep.push_back(9);
            }
            else
                espic_error("Unknown Collision Type in Process P-P Collision");
        }
    }

}

/* ------------------------------------------------------------------------- */

void Reaction::resize_threshold() 
{ 
    if(threshold.size() < static_cast<size_t>(info_size - 1)) {
        std::size_t di = info_size - 1 - threshold.size();
        while(di > 0) 
            threshold.push_back(0);  
    }
}

void Reaction::element_resize()
{
    // info_arr.resize(arr_length);
    energy.resize(arr_length);
    types.resize(info_size);
    threshold.resize(info_size);
    product_arr.resize(info_size);
}

void Reaction::find_max_coll_freq (const Real ndens)
{
    if (mr_ == 0)
        espic_error("Relative Mass not defined");
    Real e, v, nv;
    for(int i = 0; i < arr_length; ++i) {
        Real nutot = 0.;
        e = i * de_;
        v = sqrt(2.*e/mr_);
        nv = ndens * v;
        for (auto& ics :info_arr[i]){
            ics *= nv; 
            nutot += ics;
        }
        if (nutot > nu_max) { nu_max = nutot; }
    }
}

/* ------------------------------------------------------------------------- */
