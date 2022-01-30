#include "cross_section.h"

using std::vector;
using std::string;

Real kTe0;

CrossSection::CrossSection(const std::string &file)
: infile(file),
  pairs_number(0)
{
    read_input_cross_section();
    reaction_arr.reserve(num_pairs()); 
    std::cout << "Total Reaction Pair Numbers: " << pairs_number << std::endl;
    gen_reaction();

    std::cout << "Set Background Species: \n";
    std::cout << "background - (name = " << background->name;
    std::cout << ", mass = " << background->mass;
    std::cout << ", n = " << background->ndens;
    std::cout << ", vth = " << background->vth;
    std::cout << ", T = " << background->temp << ")." << std::endl;
}

CrossSection::~CrossSection()
{
    delete background;
    int num_reaction = static_cast<int>(reaction_arr.size());
    if (!reaction_arr.empty()) {
        for(int irea = 0; irea < num_reaction; ++irea) 
            delete reaction_arr[irea];
        reaction_arr.clear();
        reaction_arr.shrink_to_fit();
    }
}


    /*---------------------begin private method-----------------------*/


void CrossSection::gen_reaction() 
{ 
    std::string file;
    ReactPair spec_pair;
    Real vth = background->vth;

    for(int i = 0; i < pairs_number; ++i) {
        file = reaction_file[i];
        spec_pair = reactant_arr[i];
        reaction_arr.push_back(new Reaction(file, spec_pair, vth, i));
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
        // else if ("total" == word.at(0)) {
        //     pairs_number = atoi(word[1].c_str());
        //     word.erase(word.begin(), word.begin()+2);
        // }
        else if ("pairs" == word.at(0)) { 
            int num_re = atoi(word[1].c_str());
            word.erase(word.begin(), word.begin()+2);
            ++pairs_number;
            switch(num_re){
              case 0:
                espic_error("Insufficient Reactants");
              case 1:
                reactant_arr.push_back(std::make_pair(word[0], "default"));
                break;
              case 2:
                reactant_arr.push_back(std::make_pair(word[0], word[1]));
                break;
            }
            word.erase(word.begin(), word.begin() + num_re);
            if ("dir" == word.at(0)) {
                reaction_file.push_back(dir+word[1]);
                word.erase(word.begin(), word.begin() + 2);
            }
            else {
                    std::string cmd(word[0]);
                    espic_error(unknown_cmd_info(cmd, infile));
            }
                // word.erase(word.begin(), word.begin() + 2);
            
                // else if ("product" == word.at(0))
                // {
                //     int num_product = atoi(word[1].c_str());
                //     StringList templist;
                //     word.erase(word.begin(), word.begin()+2);
                //     if (word.size()-2 != (size_t)num_product) 
                //             espic_error("Wrong number of Products");
                //     if (num_product == 0) templist.push_back("none");
                //     else 
                //         templist.insert(templist.end(), word.begin(), word.begin()+num_product);     
                //     product_arr.push_back(std::move(templist));
                //     word.erase(word.begin(), word.begin()+num_product);
                
                // else if ("sub_cycle" == word.at(0)) {
                //     sub_cycle_arr.push_back(atoi(word[1].c_str()));
                //     word.erase(word.begin(), word.begin()+2);
                
                
            
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
    bspname = word[1];
    Real mass = (Real)atof(word[2].c_str());
    Real charge = (Real)atof(word[3].c_str());
    Real ndens = (Real)atof(word[4].c_str());
    Real temp = (Real)atof(word[5].c_str());
    background = new Background(bspname, mass, charge, ndens, temp);
}
