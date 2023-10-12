#include "octree.hpp"
#include "../simulation/simulation.hpp"

using namespace cgp;

void update_grid(ParticleArray& particles, const Octree<std::set<int>>& grid)
{
    for (int id = 0; id < particles.size(); ++id)
    {
        unsigned int octant = grid.get_octant(particles.positions[id]);

        if (octant != particles.octants[id]) {
            grid.erase_from_node(particles.octants[id], id);
            particles.octants[id] = octant;
            grid.insert_into_node(octant, id);
        }
    }
}



