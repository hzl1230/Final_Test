#include "cross_section.h"

using std::vector;
using std::string;

Real kTe0;

CrossSection::CrossSection(const std::string &file)
: infile(file) {
    read_input_cross_section();
    reaction_arr.reserve(num_pairs()); 
    get_reaction();
}

CrossSection::~CrossSection()
{
    delete background;
    for(std::size_t irea = 0; irea < reaction_arr.size(); ++irea) 
        delete reaction_arr[irea];
    reaction_arr.clear();
    reaction_arr.shrink_to_fit();
}


    /*---------------------begin private method-----------------------*/


void CrossSection::get_reaction() 
{ 
    std::string file;
    ReactPair spec_pair;
    int num_type;
    for(int i = 0; i < pairs_number; ++i) {
        StringList types;
        file = reaction_file[i];
        num_type = reaction_type_number[i];
        spec_pair = reactant_arr[i];
        types.insert(types.begin(), reaction_types.begin(), reaction_types.begin()+num_type);
        reaction_arr.emplace_back(new Reaction(file, num_type, types, spec_pair, i));
        reaction_types.erase(reaction_types.begin(), reaction_types.begin() + num_type);
    }
}

void CrossSection::read_input_cross_section()
{
    FILE *fp = fopen(infile.c_str(), "r");
    std::string dir("reaction/");
    if (NULL == fp)
    {
        std::ostringstream oss;
        oss << "Cannot read file [" << infile << "]";
        espic_error(oss.str());
    }
    std::vector<std::string> word;
    // std::string cmd;
    while (ParseLine(word, fp))
    {
        if (word.empty())
            continue;
        else if ("total" == word.at(0)) {
            pairs_number = static_cast<int>(atoi(word[1].c_str()));
        }
        else if ("pairs" == word.at(0)) { 
            int num_re = atoi(word[1].c_str());
            word.erase(word.begin(), word.begin()+2);
            switch(num_re){
              case 0:
                espic_error("Insufficient Reactants");
              case 1:
                reactant_arr.emplace_back(std::make_pair(word[0], "default"));
              case 2:
                reactant_arr.emplace_back(std::make_pair(word[0], word[1]));
            }
            word.erase(word.begin(), word.begin() + num_re);
            
            while (word.size() > 0) {
                if ("reaction" == word.at(0))
                { 
                    int nreact = atoi(std::move(word[1]).c_str());
                    reaction_type_number.emplace_back(nreact);
                    word.erase(word.begin(), word.begin() + 2);
                    for(int i = 0; i < nreact; ++i) 
                        reaction_types.emplace_back(word[i]);
                    word.erase(word.begin(), word.begin()+nreact);
                }
                else if ("dir" == word.at(0))
                {
                    reaction_file.emplace_back(dir+word[1]);
                    word.erase(word.begin(), word.begin() + 2);
                }
                else if ("product" == word.at(0))
                {
                    int num_product = atoi(word[1].c_str());
                    StringList templist;
                    word.erase(word.begin(), word.begin()+2);
                    if (word.size()-2 != (size_t)num_product) 
                            espic_error("Wrong number of Products");
                    if (num_product == 0) templist.emplace_back("none");
                    else 
                        templist.insert(templist.end(), word.begin(), word.begin()+num_product);
                        
                    product_arr.emplace_back(std::move(templist));
                    word.erase(word.begin(), word.begin()+num_product);
                }
                else
                {
                    std::string cmd(word[0]);
                    espic_error(unknown_cmd_info(cmd, infile));
                }
            }
        }
        else if ("aid_param" == word.at(0)) {
            kTe0 = (Real)atof(word[1].c_str());
        }
        else if ("background" == word.at(0)) 
            proc_background(word);
        else {
            std::string cmd(word[0]);
            espic_error(unknown_cmd_info(cmd, infile));
        }
    }
}

void CrossSection::proc_background(vector<string>& word)
{
    string cmd(word[0]);
    if (6!=word.size()) espic_error(illegal_cmd_info(cmd, infile));
    string name = word[1];
    Real mass = (Real)atof(word[2].c_str());
    Real charge = (Real)atof(word[3].c_str());
    Real ndens = (Real)atof(word[4].c_str());
    Real temp = (Real)atof(word[5].c_str());
    background = new Background(name, mass, charge, ndens, temp);
}
