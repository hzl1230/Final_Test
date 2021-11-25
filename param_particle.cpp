#include <iostream>
#include <iomanip>
#include <iterator>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>

#include "espic_math.h"
#include "parse.h"
#include "espic_info.h"
#include "species.h"
#include "ambient.h"
#include "mesh.h"
#include "param_particle.h"
#include "Inject/beam.h"
#include "Inject/inject.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ostringstream;
using std::ofstream;

/* ---------------- Begin Public Methods ---------------- */

/* Constructor */
ParamParticle::ParamParticle(const string& file, const Mesh* msh) 
  : injectdef_ptr (new InjectDef ()),
    infile(file),
    mesh(msh)
{
  init();

  cout << "Set " << num_species() << " species for particles:\n";
  for (size_type ispec = 0; ispec < num_species(); ++ispec) {
    const SpeciesDef* species = specdef_arr[ispec];
    cout << "species[" << ispec << "] - (name = " << species->name
      << " , mass = "   << species->mass
      << " , charge = " << species->charge
      << " , weight = " << species->weight << ").\n";
  }

  if (!ambientdef_arr.empty()) {
    int num_ambient = static_cast<int> (ambientdef_arr.size());

    cout << "Set " << num_ambient << " ambient environments:\n";
    for (int i = 0; i < num_ambient; i++) {
      const AmbientDef* ambient = ambientdef_arr[i];
      cout << "ambient[" << i << "] - (species = " << specdef_arr[ambient->specid]->name
        << ", n = " << ambient->ndens
        << ", T = " << ambient->temp
        << ", v = [" << ambient->vel[0]
        << ", " << ambient->vel[1]
        << ", " << ambient->vel[2] << "]).\n";
    }
  }

  if (injectdef_ptr->num_beams() > 0) {
    int num_beams = injectdef_ptr->num_beams();
    
    for (int i = 0; i < num_beams; i++) {
      const BeamDef* const& beam = injectdef_ptr->beamdef_arr[i];
      cout << "beam[" << i << "] - (species = " << specdef_arr[beam->specid]->name
        << ", n = " << beam->ndens
        << ", T = " << beam->temp
        << ", v = [" << beam->vel[0]
        << ", " << beam->vel[1]
        << ", " << beam->vel[2] << "]"
        << ", width_x = " << beam->width_x
        << ", width_y = " << beam->width_y
        << ", center = [" << beam->center[0]
        << ", " << beam->center[1]
        << ", " << beam->center[2] <<"]"
        << ", direction = [" << beam->fnorm[0]
        << ", " << beam->fnorm[1]
        << ", " << beam->fnorm[2] << "]).\n";
    }
  }

}

/* ------------------------------------------------------- */

ParamParticle::~ParamParticle()
{
  for (std::size_t iamb = 0; iamb < ambientdef_arr.size(); iamb++)
    delete ambientdef_arr[iamb];

  for (std::size_t ispec = 0; ispec < specdef_arr.size(); ispec++)
    delete specdef_arr[ispec];

  delete injectdef_ptr;
  // delete background;

  ambientdef_arr.clear();
  ambientdef_arr.shrink_to_fit();
  specdef_arr.clear();
  specdef_arr.shrink_to_fit();
}

/* ---------------- End Public Methods ---------------- */

/* ---------------- Begin Private Methods ---------------- */

/* ------------------------------------------------------- */

void ParamParticle::init()
{
  FILE *fp = fopen(infile.c_str(), "r");

  if (NULL == fp) {
    ostringstream oss;
    oss << "Cannot read file [" << infile << "]";
    espic_error(oss.str());
  }
  else 
    cout << "Read parameters to define simulatioin particle properties from [" << infile << "]" << endl;


  vector<string> word;
  while (ParseLine(word, fp)) {
    if (word.empty()) continue;     // this is a comment line

         if ("species" == word.at(0)) proc_species(word);
    else if ("ambient" == word.at(0)) proc_ambient(word);
    else if ("beam"    == word.at(0)) proc_beam(word);
    else espic_error(unknown_cmd_info(word.at(0), infile));
  }
  fclose(fp); 
}

/* ------------------------------------------------------- */

// void ParamParticle::proc_background(vector<string>& word)
// {
//   string cmd(word[0]);
//   if (6!=word.size()) espic_error(illegal_cmd_info(cmd, infile));
//   string name = word[1];
//   Real mass = (Real)atof(word[2].c_str());
//   Real charge = (Real)atof(word[3].c_str());
//   Real ndens = (Real)atof(word[4].c_str());
//   Real temp = (Real)atof(word[5].c_str());
//   background = new Background(name, mass, charge, ndens, temp);
// }

/* ------------------------------------------------------- */

void ParamParticle::proc_species(vector<string>& word)
{
  string cmd(word[0]);
  if (5 != word.size()) espic_error(illegal_cmd_info(cmd, infile));

  string name = word[1];
  Real mass = (Real)atof(word[2].c_str());
  Real charge = (Real)atof(word[3].c_str());
  Real weight = (Real)atof(word[4].c_str());
  auto search = map_spec_name_indx.find(name);
  if (search != map_spec_name_indx.end()) {
    ostringstream oss;
    oss << "species \"" << name << "\" is defined multiple times";
    espic_error(oss.str());
  }
  map_spec_name_indx[name] = num_species();
  specdef_arr.push_back(new SpeciesDef(name, mass, charge, weight));
}

/* ------------------------------------------------------- */

void ParamParticle::proc_ambient(vector<string>& word)
{
  if (0 == num_species()) {
    ostringstream oss;
    oss << "\"ambient\" command must be defined after \"spec\" in [" << infile << "]";
    espic_error(oss.str());
  }

  string cmd(word[0]);
  if (word.size() < 7) espic_error(illegal_cmd_info(cmd, infile));

  int specid = -1;
  string spec_name = word[1];
  try {
    specid = map_spec_name_indx.at(spec_name);
  }
  catch (const std::out_of_range& oor) {
    ostringstream oss;
    oss << "Unknown species \"" << spec_name << "\" given to \"ambient\" command in ["
      << infile << "]";
    espic_error(oss.str());
  }

  Real n, temp, v[3];
  // default bounds of ambient region
  Real bound_lo[3] = { mesh->xmin(), mesh->ymin(), mesh->zmin() };
  Real bound_hi[3] = { mesh->xmax(), mesh->ymax(), mesh->zmax() };

  n = (Real)atof(word[2].c_str());
  temp = (Real)atof(word[3].c_str());
  for (int c = 0; c < 3; c++) v[c] = (Real)atof(word[4+c].c_str());

  word.erase(word.begin(), word.begin()+7);
  if (word.size() > 0) {  // handle "domain" keyword
    if (2 != word.size() && 7 != word.size())
       espic_error(illegal_cmd_info("ambient", infile));

    if ("domain" != word[0])
       espic_error(illegal_cmd_info("ambient", infile));

    if (2 == word.size()) {
      if ("entire" != word[1])
       espic_error(illegal_cmd_info("ambient", infile));
    }
    else {
      for (int c = 0; c < 3; c++) {
        bound_lo[c] = (Real)atof(word[1+c].c_str());
        bound_hi[c] = (Real)atof(word[4+c].c_str());
      }
    }
  }

  ambientdef_arr.push_back(new AmbientDef(specid, n, temp, v,
                                          bound_lo, bound_hi));
  return;
}

/* ------------------------------------------------------- */

void ParamParticle::proc_beam(vector<string>& word)
{
  if (0 == num_species()) {
    ostringstream oss;
    oss << "\"beam\" command must be defined after \"spec\" in [" << infile << "]";
    espic_error(oss.str());
  }

  string cmd(word[0]);
  if (word.size() < 7) espic_error(illegal_cmd_info(cmd, infile));

  int specid = -1;
  string spec_name = word[1];
  try {
    specid = map_spec_name_indx.at(spec_name);
  }
  catch (const std::out_of_range& oor) {
    ostringstream oss;
    oss << "Unknown species \"" << spec_name << "\" given to \"beam\" command in ["
      << infile << "]";
    espic_error(oss.str());
  }

  Real n, temp, v[3];
  Real width_x = 1., width_y = 0., center[3] = {0., 0., 0.}, direction[3] = {0., 0., 0.};

  n = (Real)atof(word[2].c_str());
  temp = (Real)atof(word[3].c_str());
  for (int c = 0; c < 3; c++) v[c] = (Real)atof(word[4+c].c_str());

  word.erase(word.begin(), word.begin()+7);
  while (word.size() > 0) {  // handle other keywords
    if ("width_y" == word[0]) {
      if (word.size() > 1) {
        width_y = static_cast<Real> (atof(word[1].c_str()));
        word.erase(word.begin(), word.begin()+2);
      }
      else {
        espic_error(illegal_cmd_info(cmd, infile));
      }
    }
    else if ("width_x" == word[0]) {
      if (word.size() > 1) {
        width_x = static_cast<Real> (atof(word[1].c_str()));
        word.erase(word.begin(), word.begin()+2);
      }
      else {
        espic_error(illegal_cmd_info(cmd, infile));
      }
    }
    else if ("center" == word[0]) {
      if (word.size() > 3) {
        center[0] = static_cast<Real> (atof(word[1].c_str()));
        center[1] = static_cast<Real> (atof(word[2].c_str()));
        center[2] = static_cast<Real> (atof(word[3].c_str()));
        word.erase(word.begin(), word.begin()+4);
      }
      else {
        espic_error(illegal_cmd_info(cmd, infile));
      }
    }
    else if ("direction" == word[0]) {
      if (word.size() > 3) {
        direction[0] = static_cast<Real> (atof(word[1].c_str()));
        direction[1] = static_cast<Real> (atof(word[2].c_str()));
        direction[2] = static_cast<Real> (atof(word[3].c_str()));
        word.erase(word.begin(), word.begin()+4);
      }
      else {
        espic_error(illegal_cmd_info(cmd, infile));
      }
    }
  }   // end while (word.size() > 0)
 
  injectdef_ptr->append(new BeamDef(specid, n, temp, v, 
                                    center, direction, width_x, width_y));
}

/* ---------------- End Private Methods ---------------- */

