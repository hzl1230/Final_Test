#ifndef _REACTION_ 
#define _REACTION_ 

#include "parse.h"
#include "espic_info.h"
#include "espic_type.h"
#include "espic_math.h"
#include "collision.h"
#include "species.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <map>

typedef std::vector<std::string> StringList;
typedef std::pair<std::string,std::string> ReactPair;

enum CollisionType {ELA, EXC, ROT, VIB, ION, DIS, DISA, ISO, BACK, CEX, DR, MN, EID};
class Reaction
{
public:
    Reaction (std::string file, const ReactPair& spair, const Real vth, int id);
    ~Reaction ();

    void init_collision(const Real, const Real, const Real);

    void mcc_background_collision(const int ipart,
                                  Particle& particle, 
                                  std::vector<Species*>& species_arr, 
                                  std::vector<int>& pdep,
                                  bool e_incident);

    void mcc_particle_collision(const int ipart,
                                const int tpart,
                                Particle& incidentp, 
                                Particle& targetp,
                                std::vector<Species*>& species_arr, 
                                std::vector<int>& ipdep,
                                std::vector<int>& tpdep);  // particle(incident/fast)-particle(target/slow)

    std::vector<std::vector<int>> prodid_arr;
    bool is_background_collision;

    const int size() const { return arr_length; }
    const int isize() const { return info_size; }
    const std::vector<Real>& th() const { return threshold; }
    const ReactPair& pair() const { return spec_pair; }
    const int r_index() { return reaction_id; }
    const Real de() { return de_; }
    const std::vector<Real>& csection(int i) { return info_arr[i]; }
    const StringList& get_types() { return types; }
    const StringList& get_prod(int i) { return product_arr[i]; }
    const int sub_cycle() { return n_sub; } 
    const std::string get_file() const { return infile; }

    Real& mr() { return mr_; }
    const Real& mr() const { return mr_; }
    Real max_coll_freq() { return nu_max; }

    bool is_e_incident() { return is_e_in; }

private:
    std::string infile;
    int info_size;
    const ReactPair spec_pair;
    int reaction_id, arr_length;
    Real de_, deinv_;
    std::vector<Real> threshold;
    std::vector<Real> info;
    std::vector<std::vector<Real>> info_arr; // info_arr of every energy
    std::vector<Real> energy;                // energy_arr of total energy bin
    StringList types;
    std::vector<StringList> product_arr;
    const Real vth;
    int n_sub;
    Real mr_;
    Real nu_max;
    
    bool is_e_in;
    Collisionpair* collision;
    void resize_threshold();
    void element_resize();   
    

    void CrossSectionInterplot(const Real en, std::vector<Real>& nu)
    {
        int e_index = std::min(static_cast<int>(en*deinv_ + 0.5), arr_length - 1);
        nu.resize(info_size);
        for (int i = 0; i < info_size; ++i) 
            nu[i] = info_arr[e_index][i];
    }

    const std::vector<Real>& FetchCrossSection(const Real en)
    {
        int e_index = std::min(static_cast<int>(en*deinv_ + 0.5), arr_length - 1);
        return csection(e_index);
    }

    void find_max_coll_freq(const Real);  
};

inline Real toReal(const std::string& str) 
{
    return (Real)atof(str.c_str());
}

#endif
