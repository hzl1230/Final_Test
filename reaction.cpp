#include "reaction.h"

using std::cout;
using std::endl;

Reaction::Reaction (std::string file, int type_number, StringList type,
            const ReactPair& spair, int id) 
: infile(file),
info_size(type_number),
spec_pair(spair),
reaction_id(id),
arr_length(3),
threshold(0),
types(type)
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

    is_background_collision = true;
}    

Reaction::~Reaction()
{ }

void Reaction::resize_threshold() 
{ 
    if(threshold.size() < static_cast<size_t>(info_size - 1)) {
        std::size_t di = info_size - 1 - threshold.size();
        while(di > 0) 
            threshold.emplace_back(0);  
    }
}
