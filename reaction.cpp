#include "reaction.h"

using std::cout;
using std::endl;

Reaction::Reaction (std::string file, const ReactPair& spair, int id) 
: infile(file),
spec_pair(spair),
reaction_id(id),
threshold(0), 
mr_(0), nu_max(0)
{
    // info.reserve(info_size);
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

    is_background_collision = true;
}    

Reaction::~Reaction()
{ }

/* ------------------------------------------------------------------------- */

void Reaction::find_max_coll_freq()
{
    if (mr_ == 0)
        espic_error("Relative Mass not defined");
    Real e, v, nutot;
    for(int i = 0; i < arr_length; ++i) {
        e = i * de_;
        v = sqrt(2.*e/mr_);

        nutot = std::accumulate(info_arr[i].begin(), info_arr[i].end(), 0.);
        nutot *= v; 
        if (nutot > nu_max) { nu_max = nutot; }
    }
}

/* ------------------------------------------------------------------------- */

void Reaction::element_resize()
{
    // info_arr.resize(arr_length);
    energy.resize(arr_length);
    types.resize(info_size);
    threshold.resize(info_size);
    product_arr.resize(info_size);
}

/* ------------------------------------------------------------------------- */

void Reaction::resize_threshold() 
{ 
    if(threshold.size() < static_cast<size_t>(info_size - 1)) {
        std::size_t di = info_size - 1 - threshold.size();
        while(di > 0) 
            threshold.emplace_back(0);  
    }
}
