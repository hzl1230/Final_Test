#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <cstdio>
#include <cmath>

// #include "control.h"
#include "espic_info.h"
#include "espic_math.h"
#include "mesh.h"
#include "parse.h"

using namespace ESPIC;
using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;

Mesh::Mesh(const string& file) 
  : infile(file), 
    ndim(2), 
    nnd(-1),
    nc(-1),
    condid_field(nullptr)
{
  init();

  fbc_info[FBCType::dirichlet] = "dirichlet";
  fbc_info[FBCType::neumann]   = "neumann";
  fbc_info[FBCType::periodic]  = "periodic";
  fbc_info[FBCType::symmetric] = "symmetric";
  pbc_info[PBCType::vacuum]    = "vacuum";
  pbc_info[PBCType::reflect]   = "reflect";
  pbc_info[PBCType::periodic]  = "periodic";
//   pbc_info[symmetric] = "symmetric";

  cout << "Set simulation domain: (xmin, ymin, zmin) = (" << xmin() << ", " << ymin() << ", " << zmin() << "), ";
  cout << "(xmax, ymax, zmax) = (" << xmax() << ", " << ymax() << ", " << zmax() << ")\n";
  cout << "Set mesh: ";
  cout << "# of cells (nx, ny, nz) = (" << num_cells(0) << ", " << num_cells(1) << ", " << num_cells(2) << "), ";
  cout << "# of cells in a tile set to be (" << tile_num_cells(0)
    << ", " << tile_num_cells(1) << ", " << tile_num_cells(2) << "), ";
  cout << "cell size = (" << dx() << ", " << dy() << ", " << dz() << ").";
  cout << "\n";
  cout << "Set boundary condition for Poisson's solver:\n";
  cout << "(xmin, xmax, ymin, ymax, zmin, zmax) = (" << fbc_info.at(fbc_type(0));
  for (int k = 1; k < 6; ++k) cout << ", " << fbc_info.at(fbc_type(k));
  cout << ")\n";
  cout << "Set boundary condition for particles:\n";
  cout << "(xmin, xmax, ymin, ymax, zmin, zmax) = (" << pbc_info.at(pbc_type(0));
  for (int k = 1; k < 6; ++k) cout << ", " << pbc_info.at(pbc_type(k));
  cout << ")" << endl;

}

Mesh::~Mesh()
{
  if (condid_field != nullptr) delete [] condid_field;

  for (int i = 0; i < num_conductors(); i++) delete conductor_arr[i];
  conductor_arr.clear();
  conductor_arr.shrink_to_fit();
}

void Mesh::init()
{
  FILE *fp = fopen(infile.c_str(), "r");

  if (NULL == fp) {
    ostringstream oss;
    oss << "Cannot read file [" << infile << "]";
    espic_error(oss.str());
  }
  else 
    cout << "Read mesh control commands from [" << infile << "]" << std::endl;

  vector<string> word;
  while (ParseLine(word, fp)) {
    if (word.empty()) continue;     // this is a comment or blank line

         if ("domain"    == word.at(0)) proc_domain(word);
    else if ("num_cells" == word.at(0)) proc_num_cells(word);
    else if ("tile"      == word.at(0)) proc_tile(word);
    else if ("field_bc"  == word.at(0)) proc_field_bc(word);
    else if ("part_bc"   == word.at(0)) proc_part_bc(word);
    else if ("conductor" == word.at(0)) proc_conductor(word);
    else {
      ostringstream oss;
      oss << "Undefined command \"" << word.at(0) << "\" in [" << infile << "]";
      espic_error(oss.str());
    }
  }

  if (-1 == nnd || -1 == nc) {
    espic_error("number of mesh nodes is not defined properly");
  }

  for (int a = 0; a < 3; a++) cell_size[a] = (bound_hi[a] - bound_lo[a])/ncells[a];

  init_condid();
}

/* ------------------------------------------------------- */

void Mesh::init_condid()
{
  if (conductor_arr.empty()) return;

  condid_field = new Real[num_nodes()];
  for (Index i = 0; i < num_nodes(); i++) condid_field[i] = 0;

  Real px, py, pz;
  Index iznxy, off;
  for (int icond = 0; icond < num_conductors(); icond++) {
    Conductor* & conductor = conductor_arr[icond];
    Real condid = static_cast<Real> (conductor->get_id());

    for (Index k = 0; k < num_nodes(2); k++) {
      pz = z(0, 0, k);
      iznxy = k*num_nodes(0)*num_nodes(1);
      for (Index j = 0; j < num_nodes(1); j++) {
        py = y(0, j, 0);
        off = iznxy + j*num_nodes(0);
        for (Index i = 0; i < num_nodes(0); i++) {
          px = x(i, 0, 0);

          if (conductor->is_inside(px, py, pz)) condid_field[off + i] += condid;
        }
      }
    }

  } // end for (int icond = 0; icond < num_conductors(); icond++)

  if (dimension() == 3) {
    Index kknxy, off1, node_indx;
    
    // neigboring nodes of conductor for interior points
    for (Index k = 1; k < num_cells(2); k++) {
      iznxy = k*num_nodes(0)*num_nodes(1);

      for (Index j = 1; j < num_cells(1); j++) {
        off = iznxy + j*num_nodes(0);

        for (Index i = 1; i < num_cells(0); i++) {

          if (condid_field[off+i] > 0.5) {
            for (int kk = -1; kk < 2; kk++) {
              kknxy = (k+kk)*num_nodes(0)*num_nodes(1);
              for (int jj = -1; jj < 2; jj++) {
                off1 = kknxy + (j+jj)*num_nodes(0);
                for (int ii = -1; ii < 2; ii++) {
                  node_indx = off1 + i+ii;
                  if (condid_field[node_indx] >= 0.5) continue;
                  condid_field[node_indx] = 0.5;
                }
              }
            }
          } // end if (c_temp > 0.5)

        }
      }
    } // end for (Index k = 1; k < num_cells(2); k++)

  }
  else {
    Index off1, node_indx;

    // neigboring nodes of conductor for interior points
    for (Index j = 1; j < num_cells(1); j++) {
      off = j*num_nodes(0);

      for (Index i = 1; i < num_cells(0); i++) {

        if (condid_field[off+i] > 0.5) {
          for (int jj = -1; jj < 2; jj++) {
            off1 = (j+jj)*num_nodes(0);
            for (int ii = -1; ii < 2; ii++) {
              node_indx = off1 + i+ii;
              if (condid_field[node_indx] >= 0.5) continue;
              condid_field[node_indx] = 0.5;
            }
          }
        } // end if (c_temp > 0.5)

      }
    } // end for (Index j = 1; j < num_cells(1); j++)
    
  } // end if (dimension() == 3) else

  for (int icond = 0; icond < num_conductors(); icond++) {
    Conductor*& conductor = conductor_arr[icond];
    map_condid_arrid[conductor->get_id()] = icond;
  }
}

/* ------------------------------------------------------- */

void Mesh::proc_domain(vector<string>& word)
{
  string cmd(word[0]);
  if (7 != word.size()) espic_error(illegal_cmd_info(cmd, infile));

  auto it = word.cbegin()+1; 
  for (int a = 0; a < 3; ++a) {
    bound_lo[a] = (Real)atof((it++)->c_str());
    bound_hi[a] = (Real)atof((it++)->c_str());
  }
  if (3 != dimension()) { bound_lo[2] = 0; bound_hi[2] = 1; }
}

/* ------------------------------------------------------- */

void Mesh::proc_num_cells(vector<string>& word)
{
  string cmd(word[0]);
  if (4 != word.size()) espic_error(illegal_cmd_info(cmd, infile));

  auto it = word.cbegin()+1;
  char c[3] = {'x', 'y', 'z'};
  
  for (int a = 0; a < 3; a++) {
    ncells[a] = (Index)atoi((it++)->c_str());
    nnodes[a] = ncells[a]+1;
  }

  int n = (5 == dimension() ? 2 : dimension());
  for (int a = 0; a < n; a++) {
    if (ncells[a] < 2) { 
      ostringstream oss;
      oss << "Number of cells in the " << c[a] << " direction must be at least 2";
      espic_error(oss.str());
    }
  }
  
  if (3 != dimension()) nnodes[2] = ncells[2] = 1;  // nz = 1 for 2D or axi simulation
  nnd = nnodes[0]*nnodes[1]*nnodes[2];
  nc  = ncells[0]*ncells[1]*ncells[2];
}

/* ------------------------------------------------------- */

void Mesh::proc_tile(vector<string>& word)
{
  string cmd(word[0]);
  if (4 != word.size()) espic_error(illegal_cmd_info(cmd, infile));

  auto it = word.cbegin()+1;
  char c[3] = {'x', 'y', 'z'};
  
  for (int a = 0; a < 3; a++) {
    tncells[a] = (Index)atoi((it++)->c_str());
    tnnodes[a] = tncells[a]+1;
  }
  
  int n = (5 == dimension() ? 2 : dimension());
  for (int a = 0; a < n; a++) {
    if (tncells[a] < 1) { 
      ostringstream oss;
      oss << "Number of cells in a tile in the " << c[a] << " direction must be at least 1";
      espic_error(oss.str());
    }
  }
  
  if (3 != dimension()) tnnodes[2] = tncells[2] = 1;  // nz = 1 for 2D or axi simulation
  tnc = tncells[0]*tncells[1]*tncells[2];
  tnnd = tnnodes[0]*tnnodes[1]*tnnodes[2];
}

/* ------------------------------------------------------- */

void Mesh::proc_field_bc(vector<string>& word)
{
  FBCType b_type[6];
  Real b_val[6];

  string cmd(word[0]);
  if (word.size() < 8) espic_error(illegal_cmd_info(cmd, infile));

  /* process type */
  if ("type" != word[1]) espic_error(illegal_cmd_info(cmd, infile)); 

  word.erase(word.begin(), word.begin()+2);
  
  int i = 0;
  for (auto it = word.cbegin(); it != word.cbegin()+6; ++it) {
         if (*it == "d") b_type[i++] = FBCType::dirichlet;
    else if (*it == "n") b_type[i++] = FBCType::neumann;
    else if (*it == "p") b_type[i++] = FBCType::periodic;
    else if (*it == "s") b_type[i++] = FBCType::symmetric;
    else espic_error(illegal_cmd_info(cmd, infile));
  }
  
  // axi-symmetric, YLO must be symmetric
  if (5 == dimension() && b_type[2] != FBCType::symmetric) {
    espic_error("Field BC at YLO must be symmetric for axi-symmetric setup");
  }

  word.erase(word.begin(), word.begin()+6);

  /* process bc value */
  for (i = 0; i < 6; ++i) b_val[i] = 0.0; // default bc value - 0
  
  if (!word.empty()) {
    if ("value" != word[0]) espic_error(illegal_cmd_info(cmd, infile));
    if (7 != word.size()) espic_error(illegal_cmd_info(cmd, infile));

    for (i = 0; i < 6; ++i) b_val[i] = (Real)atof(word.at(i+1).c_str());
  }

  for (i = 0; i < 6; ++i) fbc[i] = std::make_pair(b_type[i], b_val[i]);

//  cout << "bc: type = (" << b_type[0] << ", " << b_type[1] << ", " 
//            << b_type[2] << ", " << b_type[3] << ", " 
//            << b_type[4] << ", " << b_type[5] << ")\n";
//  cout << "bc: value = (" << b_val[0] << ", " << b_val[1] << ", " 
//            << b_val[2] << ", " << b_val[3] << ", " 
//            << b_val[4] << ", " << b_val[5] << ")\n";
}

/* ------------------------------------------------------- */

void Mesh::proc_part_bc(vector<string>& word)
{
  PBCType b_type[6];
  Real b_val[6] ;

  string cmd(word[0]);
  if (word.size() < 8) espic_error(illegal_cmd_info(cmd, infile));

  /* process type */
  if ("type" != word[1]) espic_error(illegal_cmd_info(cmd, infile));

  word.erase(word.begin(), word.begin()+2);
  
  int i = 0;
  for (auto it = word.cbegin(); it != word.cbegin()+6; ++it) {
         if (*it == "v") b_type[i++] = PBCType::vacuum;
    else if (*it == "r") b_type[i++] = PBCType::reflect;
    else if (*it == "p") b_type[i++] = PBCType::periodic;
    else espic_error(illegal_cmd_info(cmd, infile));
  }
  
  // axi-symmetric, YLO must be reflection
  if (5 == dimension() && b_type[2] != PBCType::reflect) {
    espic_error("Particle BC at YLO must be reflection for axi-symmetric setup");
  }

  word.erase(word.begin(), word.begin()+6);

  /* process bc value */
  for (i = 0; i < 6; ++i) b_val[i] = 0.0; // default bc value - 0
  
  if (!word.empty()) {
    if ("value" != word[0]) espic_error(illegal_cmd_info(cmd, infile));
    if (7 != word.size()) espic_error(illegal_cmd_info(cmd, infile));

    for (i = 0; i < 6; ++i) b_val[i] = (Real)atof(word.at(i+1).c_str());
  }

  for (i = 0; i < 6; ++i) pbc[i] = std::make_pair(b_type[i], b_val[i]);

}

/* ------------------------------------------------------- */

void Mesh::proc_conductor(vector<string>& word)
{
  string cmd(word[0]);
  if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));

       if ("rectangle" == word[1]) proc_conductor_rectangle(word);
  else if ("circle"    == word[1]) proc_conductor_circle(word);
  else espic_error(illegal_cmd_info(cmd, infile));
}

/* ------------------------------------------------------- */

void Mesh::proc_conductor_rectangle(vector<string>& word)
{
  if (dimension() == 3) 
    espic_error("[conductor rectangle] works only for 2D or axi-symmetric simulations");
  string cmd(word[0]);
  ConductorRectangleDef cdef;
  bool pos_defined = false;
  bool potential_fixed = false;

  word.erase(word.begin(), word.begin()+2);

  while (!word.empty()) {
    if ("type" == word[0]) {
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      if (cdef.type != 2) { // skip if floating (must be real) already defined
        if ("virtual" == word[1]) cdef.type = 0;
        else cdef.type = 1;
      }
      word.erase(word.begin(), word.begin()+2);
    }
    else if ("epsilon" == word[0]) {
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      cdef.epsilon = static_cast<Real> (atof(word[1].c_str()));
      word.erase(word.begin(), word.begin()+2);
    }
    else if ("potential" == word[0]) {
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      if ("floating" == word[1]) {
        // An object with floating potential must be real
        // because it needs to collect charges
        cdef.type = 2;
        word.erase(word.begin(), word.begin()+2);
      }
      else if ("fixed" == word[1]) {
        if (word.size() < 3) espic_error(illegal_cmd_info(cmd, infile));
        // An object with fixed potential can be either virtual or real
//         cdef.type = 1;
        cdef.phi = static_cast<Real> (atof(word[2].c_str()));
        word.erase(word.begin(), word.begin()+3);
        potential_fixed = true;
      }
      else espic_error(illegal_cmd_info(cmd, infile));
    }
    else if ("position" == word[0]) {
      pos_defined = true;
      if (word.size() < 7) espic_error(illegal_cmd_info(cmd, infile));
      auto it = word.cbegin()+1; 
      for (int a = 0; a < 3; a++) {
        cdef.blo[a] = static_cast<Real> (atof((it++)->c_str()));
        cdef.bhi[a] = static_cast<Real> (atof((it++)->c_str()));
      }
      word.erase(word.begin(), word.begin()+7);
    }
    else if ("is_rf" == word[0]) {
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      if (word[1] == "true" && potential_fixed) cdef.rf = true;
      else cdef.rf = false;
      word.erase(word.begin(), word.begin()+2);
    }
    else espic_error(illegal_cmd_info(cmd, infile));
  }

  if (!pos_defined) {
    espic_error(illegal_cmd_info(cmd, infile));
  }

  for (int a = 0; a < 3; a++) cdef.center[a] = 0.5*(cdef.blo[a]+cdef.bhi[a]);

  conductor_arr.push_back(new ConductorRectangle(cdef));

}

/* ------------------------------------------------------- */

void Mesh::proc_conductor_circle(vector<string>& word)
{
  if (dimension() == 3) 
    espic_error("[conductor circle] works only for 2D or axi-symmetric simulations");

  string cmd(word[0]);
  ConductorCircleDef cdef;
  bool center_defined = false, radius_defined = false;

  word.erase(word.begin(), word.begin()+2);

  while (!word.empty()) {
    if ("type" == word[0]) {
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      if (cdef.type != 2) { // skip if floating (must be real) already defined
        if ("virtual" == word[1]) cdef.type = 0;
        else cdef.type = 1;
      }
      word.erase(word.begin(), word.begin()+2);
    }
    else if ("epsilon" == word[0]) {
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      cdef.epsilon = static_cast<Real> (atof(word[1].c_str()));
      word.erase(word.begin(), word.begin()+2);
    }
    else if ("potential" == word[0]) {
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      if ("floating" == word[1]) {
        // An object with floating potential must be real
        // because it needs to collect charges
        cdef.type = 2;
        word.erase(word.begin(), word.begin()+2);
      }
      else if ("fixed" == word[1]) {
        if (word.size() < 3) espic_error(illegal_cmd_info(cmd, infile));
        // An object with floating potential can be virtual or real
//         cdef.type = 1;
        cdef.phi = static_cast<Real> (atof(word[2].c_str()));
        word.erase(word.begin(), word.begin()+3);
      }
      else espic_error(illegal_cmd_info(cmd, infile));
    }
    else if ("center" == word[0]) {
      center_defined = true;
      if (word.size() < 4) espic_error(illegal_cmd_info(cmd, infile));
      auto it = word.cbegin(); 
      for (int a = 0; a < 3; a++) {
        cdef.center[a] = static_cast<Real> (atof((++it)->c_str()));
      }
      word.erase(word.begin(), word.begin()+4);
    }
    else if ("radius" == word[0]) {
      radius_defined = true;
      if (word.size() < 2) espic_error(illegal_cmd_info(cmd, infile));
      cdef.radius = static_cast<Real> (atof(word[1].c_str()));
      word.erase(word.begin(), word.begin()+2);
    }
    else espic_error(illegal_cmd_info(cmd, infile));
  }

  if ( !(center_defined && radius_defined) ) {
    espic_error(illegal_cmd_info(cmd, infile));
  }

  conductor_arr.push_back(new ConductorCircle(cdef));

}

/* ------------------------------------------------------- */

bool Mesh::is_fixed_potential(int i, int j)
{
  int condid = static_cast<int> (condid_field[j*num_nodes(0)+i]+0.0001);
  if (condid < 1) return false;
  else {
    return conductor_arr[map_condid_arrid[condid]]->is_fixed_potential();
  }
}


/* ------------------------------------------------------- */

bool Mesh::is_fixed_potential(int i, int j) const
{
  int condid = static_cast<int> (condid_field[j*num_nodes(0)+i]+0.0001);
  if (condid < 1) return false;
  else {
    return conductor_arr[map_condid_arrid.at(condid)]->is_fixed_potential();
  }
}


/* ------------------------------------------------------- */

bool Mesh::is_fixed_potential(int i, int j, int k)
{
  int condid = static_cast<int> (
    condid_field[(k*num_nodes(1)+j)*num_nodes(0)+i]+0.0001 );
  
  if (condid < 1) return false;
  else {
    return conductor_arr[map_condid_arrid[condid]]->is_fixed_potential();
  }
}

/* ------------------------------------------------------- */

bool Mesh::is_fixed_potential(int i, int j, int k) const
{
  int condid = static_cast<int> (
    condid_field[(k*num_nodes(1)+j)*num_nodes(0)+i]+0.0001 );
  
  if (condid < 1) return false;
  else {
    return conductor_arr[map_condid_arrid.at(condid)]->is_fixed_potential();
    }
}
