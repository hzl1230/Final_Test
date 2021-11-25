#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>
#include <utility>
#include <map>
#include <cassert>

// #include "utility.h"
#include "Object/conductors.h"

class Mesh {
  public:
    enum class BoundaryId { xlo, xhi, ylo, yhi, zlo, zhi};
    enum class FBCType { dirichlet, neumann, periodic, symmetric };
    enum class PBCType { vacuum, reflect, periodic };

    /* Constructors */
    /* Default constructor */
    Mesh(const std::string&);

    ~Mesh();

    /* Public methods */
    int dimension() const { return ndim; }
    Index num_nodes() const { return nnd; }
    Index num_cells() const { return nc; }
    Index num_nodes(int i) const { return nnodes[i]; }
    Index num_cells(int i) const { return ncells[i]; }
    int tile_num_cells() const { return tnc; }
    int tile_num_cells(int i) const { return tncells[i]; }
    int tile_num_nodes() const { return tnnd; }
    int tile_num_nodes(int i) const { return tnnodes[i]; }

    Real xmin() const { return bound_lo[0]; }
    Real xmax() const { return bound_hi[0]; }
    Real ymin() const { return bound_lo[1]; }
    Real ymax() const { return bound_hi[1]; }
    Real zmin() const { return bound_lo[2]; }
    Real zmax() const { return bound_hi[2]; }
//     Real dx() const { return cell_size[0]; }
//     Real dy() const { return cell_size[1]; }
//     Real dz() const { return cell_size[2]; }
//     Real dxinv() const { return 1./dx(); }
//     Real dyinv() const { return 1./dy(); }
//     Real dzinv() const { return 1./dz(); }
    Real dx(int i=0, int j=0, int k=0) const { return cell_size[0]; }
    Real dy(int i=0, int j=0, int k=0) const { return cell_size[1]; }
    Real dz(int i=0, int j=0, int k=0) const { return cell_size[2]; }
    Real dxinv(int i=0, int j=0, int k=0) const { return 1./dx(); }
    Real dyinv(int i=0, int j=0, int k=0) const { return 1./dy(); }
    Real dzinv(int i=0, int j=0, int k=0) const { return 1./dz(); }

    Real x(int i, int j, int k=0) const { return xmin() + i*cell_size[0]; }
    Real y(int i, int j, int k=0) const { return ymin() + j*cell_size[1]; }
    Real z(int i, int j, int k=0) const { return zmin() + k*cell_size[2]; }

    FBCType fbc_type(int  i) const { return fbc[i].first; }
    PBCType pbc_type(int  i) const { return pbc[i].first; }
    Real fbc_value(int i) const { return fbc[i].second; }

    int num_conductors() const {
      return static_cast<int> (conductor_arr.size());
    }

    std::vector<class Conductor*>& get_conductors() {
      return conductor_arr;
    }
    
    const std::vector<class Conductor*>& get_conductors() const {
      return conductor_arr;
    }

    class Conductor* get_conductor(int condid) {
      return conductor_arr[map_condid_arrid[condid]];
    }

    const class Conductor* get_conductor(int condid) const {
      return conductor_arr[map_condid_arrid.at(condid)];
    }

    class Conductor* get_conductor(int i, int j) {
      return get_conductor(get_condid(i, j));
    }

    const class Conductor* get_conductor(int i, int j) const {
      return get_conductor(get_condid(i, j));
    }

    int get_condid(int i, int j) const {
      return static_cast<int> (condid_field[j*num_nodes(0)+i]);
    }

    int get_condid(int i, int j, int k) const {
      return static_cast<int> (condid_field[(k*num_nodes(1)+j)*num_nodes(0)+i]);
    }

    Real* get_condid_field() { return condid_field; }

    const Real* get_condid_field() const { return condid_field; }

    // check for 2d/axi if a node's potential is fixed due to conductor's occupation
    bool is_fixed_potential(int, int); 

    bool is_fixed_potential(int, int) const;

    // check for 3d if a node's potential is fixed due to conductor's occupation
    bool is_fixed_potential(int i, int j, int k);

    bool is_fixed_potential(int i, int j, int k) const;

  private: 
    /* data member */
    std::string infile;           // name of input file for mesh definition
    int ndim;                     // dimensions (2, 3 or 5)
    Index nnd;                    // # of nodes
    Index nc;                     // # of cells
    int tnc;                      // # of cells in a tile
    int tnnd;                     // # of nodes in a tile
    Index nnodes[3];              // # of nodes in x, y and z
    Index ncells[3];              // # of cells in x, y and z
    int tncells[3];               // # of cells in x, y and z in a tile
    int tnnodes[3];               // # of nodes in x, y and z in a tile
    Real cell_size[3];
    Real bound_lo[3];
    Real bound_hi[3];             // global bounds of mesh

    std::pair<FBCType, Real> fbc[6];  // field boundary condition
                                      // (bc_type, bc_value)
    std::pair<PBCType, Real> pbc[6];  // field boundary condition

    std::map<FBCType, std::string> fbc_info;
    std::map<PBCType, std::string> pbc_info;

    Real *condid_field;
    std::vector<class Conductor*> conductor_arr;
    std::map<int, int> map_condid_arrid;

    /* Private methods */

    /* initiation */
    void init();
    void init_condid();
    void proc_domain(std::vector<std::string>&);
    void proc_num_cells(std::vector<std::string>&);
    void proc_tile(std::vector<std::string>&);
    void proc_field_bc(std::vector<std::string>&);
    void proc_part_bc(std::vector<std::string>&);

    void proc_conductor(std::vector<std::string>&);
    void proc_conductor_rectangle(std::vector<std::string>&);
    void proc_conductor_circle(std::vector<std::string>&);

};

#endif
