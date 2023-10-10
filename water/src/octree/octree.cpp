#include "octree.hpp"
#include "../simulation/simulation.hpp"

using namespace cgp;


sph_parameters_structure sph_parameters;

using s_int = std::set<int>;

void update_grid(numarray<particle_element>& particles, const Octree<std::set<int>>& grid)
{
    for (int id = 0; id < particles.size(); ++id)
    {
        particle_element& prt = particles[id];
        unsigned int code = grid.data.mortonCode(prt.p);

        if (code != prt.morton) {
            grid.erase_from_node(prt.morton, id);
            prt.morton = code;
            grid.insert_into_node(code, id);
        }
    }
}



