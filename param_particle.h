#ifndef _PARAM_PARTICLE_H
#define _PARAM_PARTICLE_H

#include <vector>
#include <string>
#include <map>

class ParamParticle {
  public:
    typedef std::size_t size_type;

    explicit ParamParticle(const std::string& file="particle.in", const class Mesh* =nullptr);
    ~ParamParticle();

    size_type num_species() const { return specdef_arr.size(); }

    std::vector<class SpeciesDef*> specdef_arr;
    std::vector<class AmbientDef*> ambientdef_arr;
    class InjectDef* injectdef_ptr;
    std::map<std::string, size_type> map_spec_name_indx;

  private:
    std::string infile;
    const class Mesh* mesh;
    void init();
    void proc_species(std::vector<std::string>&);
    void proc_ambient(std::vector<std::string>&);
    void proc_beam(std::vector<std::string>&);

};

#endif