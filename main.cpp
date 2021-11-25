#include "cross_section.h"
#include "param_particle.h"
#include "mesh.h"

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
    Mesh* mesh = new Mesh("mesh.in");
    ParamParticle* pp = new ParamParticle("particle.in", mesh);
    CrossSection* cross_section = new CrossSection("csection.in");
    for(int icps = 0; icps < cross_section->num_pairs(); ++icps)
    {
        ReactPair& pairs = cross_section->reactant_arr[icps];
        StringList& prod = cross_section->product_arr[icps];
        cout << pairs.first << " " << pairs.second << endl;
        for(const auto& pro: prod)
            cout << pro << " ";
        cout << endl;
    }

    delete cross_section;
    return 0;
}

