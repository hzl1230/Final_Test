#ifndef _REACTION_ 
#define _REACTION_ 

#include "parse.h"
#include "espic_info.h"
#include "espic_type.h"
#include "espic_math.h"
#include <algorithm>
#include <iostream>

using std::cout;
using std::endl;
typedef std::vector<std::string> StringList;
typedef std::pair<std::string,std::string> ReactPair;

inline Real toReal(const std::string& str) 
{
    return (Real)atof(str.c_str());
}

class Reaction
{
public:
    Reaction (std::string file, int type_number, StringList type,
             const ReactPair& spair, int id) 
    : infile(file),
    info_size(type_number),
    spec_pair(spair),
    types(type),
    reaction_id(id),
    arr_length(3),
    threshold(0)
    {
        info.reserve(info_size);
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
            else if("bin"==line.at(0)) {
                arr_length = atoi(line[1].c_str()); 
                info_arr.reserve(arr_length);
                energy.reserve(arr_length);
                if ("de" == line.at(2))
                    de_ = (Real)atof(line[3].c_str());
            }
            else if("threshold"==line.at(0)){
                transform(line.begin()+1, line.end(), back_inserter(threshold), toReal);
                resize_threshold();
            }
            
            else{
                energy.emplace_back((Real)atof(line[0].c_str()));
                transform(line.begin()+1, line.begin()+info_size+1, back_inserter(info), toReal);
                info_arr.emplace_back(std::move(info));
            }
        }
        cout << "Initial Particle Reaction " << reaction_id << ": " << endl;
        cout << "Reactant: [" << spec_pair.first << "," << spec_pair.second << "], "
             << "Threshold: [" ;
        for (const auto& th: threshold) 
            cout << th << " ";
        cout << "]. " ;
        cout << "de: " << de_ << ", " << "Cross Section Number: " 
             << info_size << ", Reaction Type: [";
        for (const auto& type: types)
            cout << " " << type ;
        cout << "]" << endl;
    }

    const std::vector<Real>& en_cs(Real en)
    {
        // Get cs_info through 1DLinearInterpoltationMethod
        Real deinv(1./de_);
        Real ei(en*deinv - 1) ;
        int elo = static_cast<int>(ei);
        if (elo > arr_length) {
            std::string var_name("info_arr_size");
            espic_error(out_bound_info(var_name, infile));
        }
        Real wgt[2], temp;
        // Cal node weight 
        wgt[1] = ei - elo;
        wgt[0] = 1 - wgt[1];
        
        if(!info.empty()) info.clear();
        for(int i = 0; i < info_size; ++i) {
            temp = wgt[0]*info_arr[elo][i] + wgt[1]*info_arr[elo+1][i];
            info.emplace_back(std::move(temp));
        }
        return info;
    }

    std::vector<int> productid_arr;


    const int size() const { return arr_length; }
    const int isize() const { return info_size; }
    const std::vector<Real>& th() const { return threshold; }
    const ReactPair& pair() const { return spec_pair; }
    const int r_index() { return reaction_id; }
    const Real de() { return de_; }
    const std::vector<Real>& csection(int i) { return info_arr[i]; }
    const StringList& get_types() { return types; }
    
    const std::string get_file() const { return infile; }
    
private:
    std::string infile;
    const int info_size;
    const ReactPair spec_pair;
    int reaction_id, arr_length;
    Real de_;
    std::vector<Real> threshold;
    std::vector<Real> info;
    std::vector<std::vector<Real>> info_arr; // info_arr of every energy
    std::vector<Real> energy;                // energy_arr of total energy bin
    StringList types;
    // StringList reactant_arr;
    // StringList product_arr;
    
    // std::vector<std::string> ion_product;

    void resize_threshold() 
    { 
        if(threshold.size() < static_cast<size_t>(info_size - 1)) {
            std::size_t di = info_size - 1 - threshold.size();
            while(di > 0) 
                threshold.emplace_back(0);  
        }
    }
};


#endif