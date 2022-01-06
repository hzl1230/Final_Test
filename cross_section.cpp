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
    get_reaction();

    std::cout << "Set Background Species: \n";
    std::cout << "background - (name = " << background->name;
    std::cout << ", mass = " << background->mass;
    std::cout << ", n = " << background->ndens;
    std::cout << ", T = " << background->temp << ")." << std::endl;
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
    // int num_type;
    for(int i = 0; i < pairs_number; ++i) {
        StringList types;
        file = reaction_file[i];
        spec_pair = reactant_arr[i];
        reaction_arr.emplace_back(new Reaction(file, spec_pair, i));
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
    std::string bkgspname;
    std::vector<std::string> word;
    // std::string cmd;
    while (ParseLine(word, fp))
    {
        if (word.empty())
            continue;
        else if ("background" == word.at(0)) {
            proc_background(word);
        }
        else if ("pairs" == word.at(0)) { 
            int num_re = atoi(word[1].c_str());
            word.erase(word.begin(), word.begin()+2);
            ++pairs_number;
            switch(num_re){
              case 0:
                espic_error("Insufficient Reactants");
              case 1:
                reactant_arr.emplace_back(std::make_pair(word[0], bspname));
                break;
              case 2:
                reactant_arr.emplace_back(std::make_pair(word[0], word[1]));
                break;
              case 3: 
                espic_error("Three-Body Collision not Considered");
            }
            word.erase(word.begin(), word.begin() + num_re);
            if ("dir" == word.at(0)) {
                reaction_file.emplace_back(dir+word[1]);
                // word.erase(word.begin(), word.begin() + 2);
            } else {
                std::string cmd(word[0]);
                espic_error(unknown_cmd_info(cmd, infile));
            }
        }
        else if ("aid_param" == word.at(0)) {
            kTe0 = (Real)atof(word[1].c_str());
        } else {
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
