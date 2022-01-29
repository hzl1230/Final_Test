#ifndef _CCSECTION_H
#define _CCSECTION_H
#
#include "reaction.h"
// extern Real kTe0;
class CrossSection
{
    friend class Reaction;
public:
    CrossSection(const std::string &file="csection.in");
    ~CrossSection();

    class Background {
        public:
        Background(std::string& name, Real m, Real q, Real n, Real T)
        : name(name), mass(m),
        charge(q), ndens(n), temp(T), vth(sqrt(2.*T/m))
        {}

        const std::string name;
        const Real mass;
        const Real charge;
        const Real ndens;
        const Real temp;
        const Real vth;
    };

    std::vector<class Reaction*> reaction_arr;
    std::vector<ReactPair> reactant_arr;
    
    Background* background;
    // Real kTe0;

    // const std::vector<std::string>& name_species() { return name; }
    const std::vector<std::string>& files() { return reaction_file; }
    // const std::vector<int>& num_react() const { return reaction_type_number; }
    int num_pairs() const { return pairs_number; }
    const std::string& background_name() const { return bspname; }

private:
    const std::string infile;
    int pairs_number;
    std::vector<std::string> reaction_file;
    std::string bspname;
    // std::vector<int> num_reactant;
    // std::vector<int> reaction_type_number;
    // StringList reaction_types;
    
    CrossSection();

    void gen_reaction();
    void read_input_cross_section();
    void proc_background(vector<string>& word);

};

#endif