#include "cross_section.h"
#include "param_particle.h"
#include "mesh.h"
#include "tile.h"

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
    Real dt = 0.1;
    Real curr_time = 0.;
    Real ntime = 10.;
    bool Loop = true;
    Mesh* mesh = new Mesh("mesh.in");
    ParamParticle* param_particle = new ParamParticle("particle.in", mesh);
    CrossSection* cross_section = new CrossSection("csection.in");
    Tile* tile = new Tile(mesh, param_particle, cross_section);
    cout << "MCC loop Start---------------------------------------------------------\n"; 
    while(Loop) {
        tile->ParticleCollisioninTiles(dt);
        curr_time += dt;
        if(curr_time >= ntime) 
            Loop = false;
    }
    cout << "MCC loop ends----------------------------------------------------------\n";
    delete mesh;
    delete cross_section;
    delete param_particle;
    delete tile;

    return 0;
}

