#include "octree.hpp"
#include "../simulation/simulation.hpp"

using namespace cgp;

void update_grid(numarray<particle_element>& particles, const Octree<std::set<int>>& grid)
{
    for (int id = 0; id < particles.size(); ++id)
    {
        particle_element& prt = particles[id];
        unsigned int octant = grid.get_octant(prt.p);

        if (octant != prt.octant) {
            grid.erase_from_node(prt.octant, id);
            prt.octant = octant;
            grid.insert_into_node(octant, id);
        }
    }
}



